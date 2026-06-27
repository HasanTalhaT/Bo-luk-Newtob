

# 🌌 Evren ve Güneş Sistemi Simülatörü (Universe Simulation)

Bu proje, astrofiziksel ve kozmolojik dinamikleri Python üzerinde modelleyen gelişmiş bir N-cisim (n-body) simülasyon motorudur.

**Projenin Temel Hedefi:** Klasik Newton mekaniği ile çalışan yerel gezegensel sistemlerin (Güneş Sistemi) hareketlerini, evrenin makro ölçekli dinamikleriyle (Karanlık Madde, Karanlık Enerji, Hubble Genişlemesi ve Süper Kütleli Kara Delikler) tek bir çatı altında birleştirmektir. Bu sayede, hem gezegenlerin yörünge mekaniklerini incelemek hem de evrenin büyük ölçekli yapılarının zaman içindeki evrimini analiz etmek mümkün kılınmıştır.

## ✨ Öne Çıkan Özellikler

* **Kapsamlı Fizik Motoru:** Sadece gezegenlerin kütleçekimsel etkileşimlerini değil; karanlık madde yoğunluğu, kritik yoğunluk ve karanlık enerjinin genişletici etkilerini de hesaba katan astrodinamik hesaplamalar.
* **Hibrit Ölçeklendirme:** Aynı simülasyon içerisinde hem Güneş Sistemimizdeki gezegenlerin detaylı yörünge hareketlerini hem de Samanyolu Galaksisi ve merkezindeki süper kütleli kara deliğin (Sagittarius A*) etkileşimlerini barındırır.
* **Gelişmiş Veri Toplama ve Analiz:** Simülasyon boyunca her adımda metrikleri ve yörünge verilerini kaydeder (`UniverseAnalyzer`). Veriler, Numpy dizilerini destekleyen özel bir JSON kodlayıcı ile güvenle dışa aktarılır.
* **Görselleştirme Altyapısı:** Matplotlib/grafik kütüphaneleri entegrasyonuna hazır, iz (trail) uzunluğu ve renk haritası gibi detayların ayarlanabildiği bir görselleştirme konfigürasyonuna sahiptir.

## 🏗️ Proje Mimarisi

Simülasyon, modüler ve nesne yönelimli (OOP) bir yaklaşımla tasarlanmıştır. Temel bileşenler şunlardır:

* `simulation.py`: Simülasyonun kalbi olan `UniverseSimulation` sınıfını ve konfigürasyon yöneticisini içerir.
* `particles.py`: Kütle, pozisyon ve hız vektörlerine sahip standart gök cisimlerini (Güneş, Dünya, Jüpiter vb.) tanımlar.
* `galaxy.py` & `blackhole.py`: Sgr A* gibi devasa yapıların ve galaktik özelliklerin modellendiği modüller.
* `constants.py`: Güneş kütlesi, gezegen kütleleri ve kozmolojik sabitleri (Hubble, Dark Matter vb.) barındıran sabitler kütüphanesi.
* `universe_analyzer.py`: Çıktıların bilimsel analizini ve yapı oluşumlarının incelenmesini sağlayan veri analiz modülü.

## 🚀 Hızlı Başlangıç

### Gereksinimler

Projeyi çalıştırmak için aşağıdaki Python kütüphanelerinin yüklü olması gerekmektedir:

```bash
pip install numpy tqdm

```

### Kullanım

Simülasyonu başlatmak için ana dosyayı çalıştırmanız yeterlidir. Varsayılan olarak simülasyon, Güneş Sistemi gök cisimleri ile 10 yıllık bir süreci (günlük zaman adımlarıyla) hesaplayacaktır:

```bash
python main.py

```

Simülasyon tamamlandığında, analiz sonuçları ekrana yazdırılacak ve yörünge metrikleri `simulation_results` klasörüne kaydedilecektir.

---

🚧 Geliştirme Durumu (Work in Progress)
Önemli Not: Bu proje şu anda aktif geliştirme aşamasındadır ve henüz tam stabil olarak çalışmamaktadır.
Gezegensel ölçek ile galaktik ölçekteki kütleçekimsel kuvvetlerin (Karanlık Enerji, Sagittarius A* gibi) aynı anda simüle edilmesi oldukça karmaşık bir matematiksel altyapı gerektirmektedir. Bu nedenle simülasyonu çalıştırdığınızda bazı fiziksel sapmalar, performans darboğazları veya beklenmeyen hatalar (bug'lar) ile karşılaşabilirsiniz.
Önemli Not: Bu oldukça eski bir projedir (legacy project) ve aktif olarak bakımı yapılmamaktadır. Yazıldığı dönemin standartlarını ve kütüphane sürümlerini yansıtmaktadır.

Bununla birlikte, gezegensel ölçek ile galaktik ölçekteki kütleçekimsel kuvvetlerin (Karanlık Enerji, Sagittarius A* vb.) aynı anda simüle edilmesi gibi zorlu bir hedefi olduğundan, proje tam stabil olarak çalışmamaktadır. Simülasyonu çalıştırdığınızda bazı fiziksel sapmalar, performans darboğazları veya beklenmeyen hatalar (bug'lar) ile karşılaşabilirsiniz.

Projenin temel amacı konsepti test etmek ve altyapıyı kurmaktır. Sorunları tespit eden veya algoritmaları iyileştirmek isteyen herkesin Pull Request (PR) katkılarına açığım!
