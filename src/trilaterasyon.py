"""Adım 5: Üç alıcıdan RSSI mesafeleriyle konum (trilaterasyon) ve harita.

Konum, sum_i (||p - p_i|| - d_i)^2 toplamını küçülten p'dir (scipy.optimize.least_squares).
Koordinatlar yerel düzlemde metre (x doğu, y kuzey). Harita için keyfi bir başlangıç noktasına
göre enlem ve boylama çevrilir; bu nokta gerçek bir arama alanını temsil etmez.

Çıktılar (results/):
  trilaterasyon.csv     geometri ve süzme durumuna göre konum hatası
  trilaterasyon.png     kestirilen konumların dağılımı
  harita.html           alıcılar, gerçek konum, kestirilen konum ve hata dairesi (folium)
"""
import folium
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import pandas as pd
from scipy.optimize import least_squares

from kalman import kalman
from ortak import SONUC, TOHUM, mesafe_kestir, rssi_uret

N_YOL = 2.7
SIGMA = 6.0
PAKET = 20          # süzülmüş durumda alıcı başına paket sayısı
KOSU = 2000
HEDEF = np.array([350.0, 250.0])
GEOMETRILER = {
    "ucgen": np.array([[0.0, 0.0], [800.0, 0.0], [400.0, 693.0]]),
    "neredeyse_dogrusal": np.array([[0.0, 0.0], [400.0, 40.0], [800.0, 0.0]]),
    "dogrusal": np.array([[0.0, 0.0], [400.0, 0.0], [800.0, 0.0]]),
}
BASLIKLAR = {"ucgen": "triangle", "neredeyse_dogrusal": "nearly collinear (middle 40 m off)",
             "dogrusal": "collinear receivers"}
BASLANGIC_ENLEM, BASLANGIC_BOYLAM = 39.0, 35.0  # keyfi harita başlangıcı (derece)


def konum_bul(alicilar, d, rng):
    """Mesafelerden konum (m). Birden çok başlangıçtan çözüp en düşük maliyetliyi seçer.

    Alıcılar aynı doğru üzerindeyse doğrunun iki yanındaki ayna noktaları tam aynı maliyeti verir
    ve veri hangisinin doğru olduğunu söyleyemez. Eşit maliyetli ayrı çözümler arasından rastgele
    seçilir; sabit bir sıra seçmek kestiriciye gerçek konumu gizlice söylemek olurdu.
    """
    merkez = alicilar.mean(axis=0)
    kalinti = lambda p: np.linalg.norm(alicilar - p, axis=1) - d
    cozumler = [least_squares(kalinti, merkez + np.array(k, dtype=float))
                for k in ([0, 0], [300, 300], [300, -300], [-300, 300], [-300, -300])]
    en_dusuk = min(s.cost for s in cozumler)
    adaylar = []
    for s in cozumler:
        ayri = all(np.linalg.norm(s.x - a) > 1.0 for a in adaylar)
        if s.cost <= en_dusuk * (1 + 1e-6) + 1e-9 and ayri:
            adaylar.append(s.x)
    return adaylar[rng.integers(len(adaylar))]


def mesafeler(alicilar, rng, paket):
    """Her alıcı için RSSI'dan mesafe (m). paket > 1 ise RSSI önce Kalman (Q = 0) ile süzülür."""
    gercek = np.linalg.norm(alicilar - HEDEF, axis=1)
    z = rssi_uret(np.repeat(gercek[:, None], paket, axis=1), N_YOL, SIGMA, rng)
    rssi = z[:, 0] if paket == 1 else kalman(z, SIGMA**2)[0][:, -1]
    return mesafe_kestir(rssi, N_YOL)


def benzet(ad, alicilar, paket, rng):
    kestirimler = np.array([konum_bul(alicilar, mesafeler(alicilar, rng, paket), rng) for _ in range(KOSU)])
    hata = np.linalg.norm(kestirimler - HEDEF, axis=1)
    ozet = {
        "geometri": ad, "paket": paket,
        "medyan_m": float(np.median(hata)), "p90_m": float(np.percentile(hata, 90)),
        "rmse_m": float(np.sqrt(np.mean(hata**2))),
        # Hedef y > 0 tarafında. Doğrusal geometride y < 0 ayna çözümüdür.
        "y_eksi_tarafi_orani": float(np.mean(kestirimler[:, 1] < 0)),
    }
    return ozet, kestirimler


def enlem_boylam(xy):
    """Yerel metre -> (enlem, boylam), küçük alan için düz yaklaşım."""
    xy = np.atleast_2d(xy)
    enlem = BASLANGIC_ENLEM + xy[:, 1] / 111_320.0
    boylam = BASLANGIC_BOYLAM + xy[:, 0] / (111_320.0 * np.cos(np.radians(BASLANGIC_ENLEM)))
    return np.column_stack([enlem, boylam])


def harita(ornekler, ozet):
    m = folium.Map(location=enlem_boylam(HEDEF)[0].tolist(), zoom_start=15, control_scale=True)
    renk = {"ucgen": "blue", "neredeyse_dogrusal": "orange", "dogrusal": "red"}
    for ad, (alicilar, d, kestirim) in ornekler.items():
        p90 = ozet[(ozet["geometri"] == ad) & (ozet["paket"] == PAKET)]["p90_m"].iloc[0]
        katman = folium.FeatureGroup(name=f"{ad} geometry ({PAKET} packets, Kalman)", show=(ad == "ucgen"))
        for i, (a, di) in enumerate(zip(enlem_boylam(alicilar), d)):
            folium.Marker(a.tolist(), tooltip=f"receiver {i + 1}, estimated distance {di:.0f} m",
                          icon=folium.Icon(color="gray", icon="signal")).add_to(katman)
            folium.Circle(a.tolist(), radius=float(di), color="gray", weight=1, fill=False).add_to(katman)
        tahmin = enlem_boylam(kestirim)[0].tolist()
        hata = float(np.linalg.norm(kestirim - HEDEF))
        folium.Marker(tahmin, tooltip=f"estimated position, error {hata:.0f} m",
                      icon=folium.Icon(color=renk[ad], icon="screenshot")).add_to(katman)
        folium.Circle(tahmin, radius=p90, color=renk[ad], fill=True, fill_opacity=0.1,
                      tooltip=f"90% error circle, r = {p90:.0f} m (simulation)").add_to(katman)
        katman.add_to(m)
    folium.Marker(enlem_boylam(HEDEF)[0].tolist(), tooltip="true position",
                  icon=folium.Icon(color="green", icon="user")).add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)
    m.save(str(SONUC / "harita.html"))


def main():
    SONUC.mkdir(exist_ok=True)
    rng = np.random.default_rng(TOHUM)

    for ad, alicilar in GEOMETRILER.items():
        hatalar = []
        for _ in range(20):
            p = konum_bul(alicilar, np.linalg.norm(alicilar - HEDEF, axis=1), rng)
            hatalar.append(np.linalg.norm(p - HEDEF))
        hatalar = np.array(hatalar)
        print(f"gürültüsüz, {ad}, 20 çözüm: en küçük hata {hatalar.min():.2e} m, "
              f"ayna noktasına ({HEDEF[0]:g}, {-HEDEF[1]:g}) giden {int(np.sum(hatalar > 1))}/20")

    satirlar, bulutlar = [], {}
    for ad, alicilar in GEOMETRILER.items():
        for paket in (1, PAKET):
            ozet, kestirimler = benzet(ad, alicilar, paket, rng)
            satirlar.append(ozet)
            bulutlar[(ad, paket)] = kestirimler
    ozet = pd.DataFrame(satirlar)
    ozet.to_csv(SONUC / "trilaterasyon.csv", index=False)
    print(f"hedef {HEDEF}, n = {N_YOL}, sigma = {SIGMA:g} dB, {KOSU} koşu")
    print(ozet.round(3).to_string(index=False))

    fig, eksenler = plt.subplots(1, 3, figsize=(17, 5.5))
    for ax, (ad, alicilar) in zip(eksenler, GEOMETRILER.items()):
        for paket, r in ((1, "0.7"), (PAKET, "#1f77b4")):
            b = bulutlar[(ad, paket)]
            ax.plot(b[:, 0], b[:, 1], ".", ms=2, color=r, alpha=0.5,
                    label="1 packet" if paket == 1 else f"{PAKET} packets, Kalman")
        ax.plot(*alicilar.T, "k^", ms=10, label="receivers")
        ax.plot(*HEDEF, "g*", ms=15, label="true position")
        ax.set_title(BASLIKLAR[ad])
        ax.set_xlim(-600, 1400)
        ax.set_ylim(-900, 1100)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.grid(alpha=0.3)
        # Bulut noktaları çok küçük; göstergede okunur boyutta temsilciler kullanılır
        temsilci = [Line2D([], [], ls="", marker="o", ms=5, color="0.7", label="1 packet"),
                    Line2D([], [], ls="", marker="o", ms=5, color="#1f77b4", label=f"{PAKET} packets, Kalman"),
                    Line2D([], [], ls="", marker="^", ms=8, color="k", label="receivers"),
                    Line2D([], [], ls="", marker="*", ms=12, color="g", label="true position")]
        ax.legend(handles=temsilci, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(SONUC / "trilaterasyon.png", dpi=150)
    plt.close(fig)

    # Harita için her geometriden tek bir süzülmüş örnek
    ornekler = {}
    for ad, alicilar in GEOMETRILER.items():
        d = mesafeler(alicilar, rng, PAKET)
        ornekler[ad] = (alicilar, d, konum_bul(alicilar, d, rng))
        print(f"harita örneği, {ad}: mesafeler {d.round(1)}, kestirim {ornekler[ad][2].round(1)}")
    harita(ornekler, ozet)
    print("yazıldı: results/trilaterasyon.csv, trilaterasyon.png, harita.html")


if __name__ == "__main__":
    main()
