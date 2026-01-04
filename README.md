# PSK Transceiver with Physical Layer Security (Artificial Noise)

This repository contains a full PSK transceiver chain implemented in MATLAB, focusing on synchronization and **Physical Layer Security (PLS)**. It demonstrates how artificial noise can protect data from unauthorized eavesdroppers.

---

## English Description

### System Architecture
The system utilizes a **Zadoff-Chu sequence** for frame synchronization and **Meyr-Oeder algorithms** for blind timing recovery. It includes a complete transmitter/receiver chain with Root-Raised Cosine (RRC) filtering and CFO correction.



### Key Features & Contributions
- **Artificial Noise (AN):** Implementation of Gaussian snd Phase-shift oise injection also XOR-based cryptography technique to degrade unauthorized receivers while maintaining performance for authorized users.
- **Robust Synchronization:** Integration of Zadoff-Chu correlation for precise packet detection and Farrow fractional interpolation for timing correction.
- **Performance Evaluation:** Comparative Bit Error Rate (BER) analysis between Authorized (Keyed) and Unauthorized (No Key) receivers.

### Academic Credit & Integrity
- The core synchronization files starting with `THAL_` are provided by the **Istanbul Technical University Wireless Communication Research Laboratory (THAL)**.
- The system integration, security layer (AN/XOR), and PSK object-oriented logic were developed by me as part of my academic work.

---

## Fiziksel Katman Güvenliği (Yapay Gürültü) Özellikli PSK Alıcı-Verici

### Sistem Mimarisi
Sistem, çerçeve senkronizasyonu için **Zadoff-Chu dizilerini** ve kör zamanlama geri kazanımı için **Meyr-Oeder algoritmalarını** kullanmaktadır. SRRC filtreleme ve CFO (Frekans Kayması) düzeltme içeren tam bir alıcı-verici zinciri sunar.

### Temel Özellikler ve Katkılar
- **Yapay Gürültü (Artificial Noise - AN):** Yetkisiz alıcıların başarımını düşürmek için Gauss ve Faz Kaydırmalı gürültü enjeksiyonu ayrıca XOR tabanlı kriptolama yöntemleri uygulanmıştır.
- **Güçlü Senkronizasyon:** Hassas paket tespiti için Zadoff-Chu korelasyonu ve zamanlama düzeltme için Farrow fraksiyonel interpolasyonu entegre edilmiştir.
- **Performans Değerlendirmesi:** Yetkili (Anahtarlı) ve yetkisiz (Anahtarsız) alıcılar arasında karşılaştırmalı Bit Hata Oranı (BER) analizi sunulmaktadır.

### Akademik Etik ve Teşekkür
- `THAL_` ön eki ile başlayan dosyalar **İstanbul Teknik Üniversitesi Telsiz Haberleşme Araştırma Laboratuvarı (THAL)** tarafından sağlanmıştır.
- Sistem entegrasyonu, güvenlik katmanının kurgulanması (AN ve XOR mekanizmaları) ve PSK nesne yönelimli mantığı akademik çalışmalarım kapsamında tarafımdan geliştirilmiştir.

---

## How to Run / Nasıl Çalıştırılır
1. Clone the repository / Depoyu klonlayın.
2. Add the `src` folder to your MATLAB path / `src` klasörünü MATLAB yoluna ekleyin.
3. Run `tests/performance_analysis.m` to see the results / Sonuçları görmek için testi çalıştırın.

---

## Requirements / Gereksinimler

To run this simulation, the following MATLAB Toolboxes are required:
/ Bu simülasyonu çalıştırmak için aşağıdaki MATLAB Toolbox'larının yüklü olması gerekir:

- **Communications Toolbox** 
- **Signal Processing Toolbox** 
