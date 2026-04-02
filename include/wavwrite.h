#ifndef WAV_WRITER_H
#define WAV_WRITER_H

#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <iostream>

struct WavHeader {
    char riffId[4] = {'R', 'I', 'F', 'F'};
    uint32_t fileSize;
    char waveId[4] = {'W', 'A', 'V', 'E'};
    char fmtId[4] = {'f', 'm', 't', ' '};
    uint32_t fmtSize = 16;
    uint16_t audioFormat = 1; // PCM
    uint16_t numChannels = 1; // Mono
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample = 16;
    char dataId[4] = {'d', 'a', 't', 'a'};
    uint32_t dataSize;
};

class WavWriter {
public:
    // data: 音圧や流量の時系列データ
    // dt: シミュレーションの時間刻み (s)
    // filename: 保存ファイル名
    static void save(const std::vector<double>& data, double dt, const std::string& filename) {
        if (data.empty()) return;

        // サンプリング周波数
        uint32_t sampleRate = static_cast<uint32_t>(std::round(1.0 / dt));
        
        // データの最大絶対値を見つけて正規化
        double maxVal = 0.0;
        for (double v : data) {
            maxVal = std::max(maxVal, std::abs(v));
        }

        if (maxVal < 1e-12) maxVal = 1.0; // ゼロ除算防止

        std::vector<int16_t> pcmData;
        pcmData.reserve(data.size());

        // double -> int16_t (-32767 ~ 32767) に変換
        
        double scale = 32760.0 / maxVal; 
        
        for (double v : data) {
            pcmData.push_back(static_cast<int16_t>(v * scale));
        }

        // ヘッダー作成
        WavHeader header;
        header.sampleRate = sampleRate;
        header.bitsPerSample = 16;
        header.numChannels = 1;
        header.blockAlign = header.numChannels * (header.bitsPerSample / 8);
        header.byteRate = header.sampleRate * header.blockAlign;
        header.dataSize = static_cast<uint32_t>(pcmData.size() * sizeof(int16_t));
        header.fileSize = 36 + header.dataSize;

        // バイナリ書き出し
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            std::cerr << "Cannot create wav file: " << filename << std::endl;
            return;
        }

        file.write(reinterpret_cast<const char*>(&header), sizeof(WavHeader));
        file.write(reinterpret_cast<const char*>(pcmData.data()), pcmData.size() * sizeof(int16_t));
        
        std::cout << "[WavWriter] Saved " << filename << " (Fs=" << sampleRate << "Hz, MaxVal=" << maxVal << ")" << std::endl;
    }
};

#endif