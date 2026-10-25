// Included only by the isolated diagnostic copy, after envelopeEdge().
// Contour, yOnEdge(), cutTolerance and standard headers come from that TU.
// Nothing here runs until the original interval-intersection check throws.

void sectionDiagnosticString(std::ostream& out, const std::string& value) {
    const char* hex = "0123456789abcdef";
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 0x20) out << "\\u00" << hex[c >> 4] << hex[c & 15];
        else out << static_cast<char>(c);
    }
    out << '"';
}

void sectionDiagnosticNumber(std::ostream& out, double value) {
    if (std::isfinite(value)) out << value;
    else out << "null";
}

void sectionDiagnosticPoint(std::ostream& out, const CutPoint& p) {
    out << "{\"y\":";
    sectionDiagnosticNumber(out, p.y);
    out << ",\"z\":";
    sectionDiagnosticNumber(out, p.z);
    out << '}';
}

void sectionDiagnosticContour(std::ostream& out, const Contour& contour,
                              double z, bool takeMaximum) {
    out << "{\"take_maximum_y\":" << (takeMaximum ? "true" : "false")
        << ",\"cut_nodes\":[";
    for (size_t i = 0; i < contour.nodes.size(); ++i) {
        if (i) out << ',';
        sectionDiagnosticPoint(out, contour.nodes[i]);
    }
    out << "],\"edges\":[";
    for (size_t i = 0; i < contour.edges.size(); ++i) {
        if (i) out << ',';
        out << '[' << contour.edges[i].a << ',' << contour.edges[i].b << ']';
    }
    out << "],\"intersections\":[";
    std::vector<double> distinct;
    bool first = true;
    for (int i = 0; i < static_cast<int>(contour.edges.size()); ++i) {
        const auto edge = contour.edges[i];
        const auto a = contour.nodes[edge.a], b = contour.nodes[edge.b];
        // Repeat the original strict midpoint test and dedup order exactly.
        if (!(z > std::min(a.z, b.z) && z < std::max(a.z, b.z))) continue;
        const double y = yOnEdge(contour, i, z);
        bool duplicate = false;
        for (double other : distinct)
            if (std::abs(y - other) <= cutTolerance) duplicate = true;
        if (!duplicate) distinct.push_back(y);
        if (!first) out << ',';
        first = false;
        out << "{\"edge_id\":" << i << ",\"cut_node_ids\":[" << edge.a << ',' << edge.b
            << "],\"a\":";
        sectionDiagnosticPoint(out, a);
        out << ",\"b\":";
        sectionDiagnosticPoint(out, b);
        out << ",\"y\":";
        sectionDiagnosticNumber(out, y);
        out << ",\"duplicate\":" << (duplicate ? "true" : "false") << '}';
    }
    out << "],\"distinct_count\":" << distinct.size() << ",\"distinct_y\":[";
    for (size_t i = 0; i < distinct.size(); ++i) {
        if (i) out << ',';
        sectionDiagnosticNumber(out, distinct[i]);
    }
    out << "]}";
}

void reportChannelSectionFailure(const Contour& left, const Contour& right,
                                 const char* failedSide, int section, double x,
                                 double z0, double z1, double mid, bool maxLeft,
                                 const char* reason) {
    std::ostringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);
    out << "{\"schema_version\":1,\"failed_side\":";
    sectionDiagnosticString(out, failedSide);
    out << ",\"reason\":";
    sectionDiagnosticString(out, reason);
    out << ",\"section\":" << section << ",\"x_mm\":";
    sectionDiagnosticNumber(out, x);
    out << ",\"z_interval_mm\":[";
    sectionDiagnosticNumber(out, z0);
    out << ',';
    sectionDiagnosticNumber(out, z1);
    out << "],\"z_mid_mm\":";
    sectionDiagnosticNumber(out, mid);
    out << ",\"dedup_tolerance_mm\":" << cutTolerance
        << ",\"node_id_space\":\"cut contour indices, not VTU point IDs\",\"left\":";
    sectionDiagnosticContour(out, left, mid, maxLeft);
    out << ",\"right\":";
    sectionDiagnosticContour(out, right, mid, !maxLeft);
    out << '}';
    std::cerr << "CHANNEL_SECTION_DIAGNOSTIC_JSON " << out.str() << '\n';
}
