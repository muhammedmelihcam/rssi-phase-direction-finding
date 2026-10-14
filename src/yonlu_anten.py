"""Adım 2: Yönlü anteni döndürerek RSSI ile yön kestirimi.

Anten kazancı dB cinsinden Gauss biçiminde: G(theta) = -12 (theta / B)^2, B 3 dB ışın genişliği (derece).
Gerçek antenlerde arka ve yan loblar sonsuz zayıf değildir; kazanç TABAN_DB altına inmez.

Üç kestirim yöntemi:
  tepe    en yüksek RSSI ölçülen açı (çözünürlük adım kadar)
  parabol en yüksek nokta ve iki komşusundan geçen parabolün tepesi. G dB cinsinden açının
          karesiyle değiştiği için parabol dB'de tam modeldir; doğrusal güçte Gauss uydurmakla aynıdır.
  gauss   ana lobdaki daha çok noktaya en küçük karelerle aynı parabol (doğrusal güçte Gauss)

Çıktılar (results/):
  yon_hata.csv           her (B, sigma, yöntem) için hata ölçütleri
  yon_hata.png           hata - gürültü ve hata - ışın genişliği
  yon_ornek_tarama.png   tek bir dönüşte ölçülen RSSI ve kestirimler
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ortak import SONUC, TOHUM, aci_farki

ADIM_DEG = 10.0                                  # anten dönme adımı
ACILAR = np.arange(0.0, 360.0, ADIM_DEG)         # 36 ölçüm açısı
TABAN_DB = -30.0                                 # yan ve arka lob tabanı (ön-arka oranı 30 dB)
P_ALINAN_DBM = -80.0                             # ana lob tepesinde alınan güç; yön kestirimini etkilemez
TEKRAR = 3000                                    # her durum için rastgele gerçek yön sayısı
ISIN_GENISLIKLERI = [20.0, 30.0, 45.0, 60.0, 90.0]
SIGMALAR = [0.0, 1.0, 2.0, 4.0, 6.0]
YONTEMLER = ["tepe", "parabol", "gauss"]


def kazanc_db(theta_deg, B_deg):
    """Anten kazancı (dB). theta, anten ekseniyle hedef arasındaki açı (derece)."""
    theta = aci_farki(theta_deg, 0.0)
    return np.maximum(-12.0 * (theta / B_deg) ** 2, TABAN_DB)


def tarama_uret(gercek_deg, B_deg, sigma_db, rng):
    """Anten ACILAR boyunca dönerken ölçülen RSSI (dBm). gercek_deg (N,) -> sonuç (N, 36)."""
    gercek = np.asarray(gercek_deg, dtype=float)[:, None]
    temiz = P_ALINAN_DBM + kazanc_db(ACILAR[None, :] - gercek, B_deg)
    return temiz + rng.normal(0.0, sigma_db, temiz.shape)


def kestir_tepe(rssi):
    """En yüksek RSSI'nın açısı (derece)."""
    return ACILAR[np.argmax(rssi, axis=1)]


def kestir_parabol(rssi):
    """Tepe noktası ve iki komşusundan geçen parabolün tepesi (derece).

    y(-1), y(0), y(+1) için tepe kayması: delta = (y(-1) - y(+1)) / (2 (y(-1) - 2 y(0) + y(+1))),
    birimi adım. Açılar daireseldir, komşular dizinin başından ve sonundan sarılır.
    """
    satir = np.arange(len(rssi))
    i = np.argmax(rssi, axis=1)
    y0 = rssi[satir, i]
    ysol = rssi[satir, (i - 1) % rssi.shape[1]]
    ysag = rssi[satir, (i + 1) % rssi.shape[1]]
    payda = ysol - 2 * y0 + ysag
    # payda < 0 tepe için doğal durumdur; düz ya da ters eğride tepe açısında kal
    gecerli = payda < 0
    delta = np.where(gecerli, (ysol - ysag) / (2 * np.where(gecerli, payda, -1.0)), 0.0)
    delta = np.clip(delta, -0.5, 0.5)  # tepe en yüksek ölçümün yarım adımı içinde olmalı
    return (ACILAR[i] + delta * ADIM_DEG) % 360.0


def kestir_gauss(rssi, B_deg):
    """Ana lobdaki 2K + 1 noktaya en küçük karelerle parabol (dB'de) uydurur, tepesini döner.

    K = round(B / adım): tepe ölçümünün iki yanında yaklaşık bir ışın genişliği kadar nokta.
    Bu noktalar -12 dB civarında kalır, TABAN_DB'ye değmez; model orada hâlâ tam paraboldür.
    Doğrusal güçte Gauss uydurmakla eşdeğerdir, ama daha çok nokta kullandığı için gürültüyü ortalar.
    """
    K = max(1, int(round(B_deg / ADIM_DEG)))
    u = np.arange(-K, K + 1)                        # tepe ölçümüne göre adım cinsinden konum
    A = np.column_stack([np.ones_like(u), u, u**2])  # y = a + b u + c u^2
    i = np.argmax(rssi, axis=1)
    sutun = (i[:, None] + u[None, :]) % rssi.shape[1]
    Y = np.take_along_axis(rssi, sutun, axis=1)
    _, b, c = np.linalg.pinv(A) @ Y.T               # bütün satırlar tek çarpımda
    gecerli = c < 0
    delta = np.where(gecerli, -b / (2 * np.where(gecerli, c, -1.0)), 0.0)
    delta = np.clip(delta, -K, K)
    return (ACILAR[i] + delta * ADIM_DEG) % 360.0


KESTIRICILER = {
    "tepe": lambda rssi, B: kestir_tepe(rssi),
    "parabol": lambda rssi, B: kestir_parabol(rssi),
    "gauss": kestir_gauss,
}


def hata_olc(B_deg, sigma_db, rng, tekrar=TEKRAR):
    gercek = rng.uniform(0.0, 360.0, tekrar)
    rssi = tarama_uret(gercek, B_deg, sigma_db, rng)
    satirlar = []
    for ad, f in KESTIRICILER.items():
        hata = np.abs(aci_farki(f(rssi, B_deg), gercek))
        satirlar.append({
            "B_deg": B_deg, "sigma_db": sigma_db, "yontem": ad,
            "rmse_deg": float(np.sqrt(np.mean(hata**2))),
            "medyan_deg": float(np.median(hata)),
            "p90_deg": float(np.percentile(hata, 90)),
            "max_deg": float(hata.max()),
        })
    return satirlar


def grafikler(df, rng):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    renk = {"tepe": "#d62728", "parabol": "#1f77b4", "gauss": "#2ca02c"}
    etiket = {"tepe": "max RSSI angle", "parabol": "parabola (3 points)",
              "gauss": "Gaussian fit (main lobe, least squares)"}
    for y in YONTEMLER:
        k = df[(df["B_deg"] == 30.0) & (df["yontem"] == y)]
        a1.plot(k["sigma_db"], k["rmse_deg"], "o-", color=renk[y], label=etiket[y])
        k = df[(df["sigma_db"] == 2.0) & (df["yontem"] == y)]
        a2.plot(k["B_deg"], k["rmse_deg"], "o-", color=renk[y], label=etiket[y])
    a1.axhline(ADIM_DEG / np.sqrt(12), color="0.5", ls=":", label="step / √12 (quantisation)")
    a1.set_xlabel("RSSI noise σ [dB]")
    a1.set_title(f"B = 30°, step = {ADIM_DEG:g}°")
    a2.set_xlabel("3 dB beamwidth B [deg]")
    a2.set_title(f"σ = 2 dB, step = {ADIM_DEG:g}°")
    for ax in (a1, a2):
        ax.set_ylabel("bearing RMSE [deg]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(SONUC / "yon_hata.png", dpi=150)
    plt.close(fig)

    # Tek bir dönüş örneği
    gercek, B, sigma = 127.0, 30.0, 2.0
    rssi = tarama_uret([gercek], B, sigma, rng)
    ince = np.linspace(0, 360, 721)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(ince, P_ALINAN_DBM + kazanc_db(ince - gercek, B), "k-", lw=1, label="noise-free pattern")
    ax.plot(ACILAR, rssi[0], "o", color="#555555", label=f"measured, σ = {sigma:g} dB")
    ax.axvline(gercek, color="k", ls=":", label=f"true {gercek:g}°")
    for y, cizgi in zip(YONTEMLER, ["-", "--", "-."]):
        kestirim = KESTIRICILER[y](rssi, B)[0]
        ax.axvline(kestirim, color=renk[y], ls=cizgi, label=f"{etiket[y]}: {kestirim:.1f}°")
    ax.set_xlabel("antenna pointing angle [deg]")
    ax.set_ylabel("RSSI [dBm]")
    ax.set_title(f"One rotation, B = {B:g}°")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(SONUC / "yon_ornek_tarama.png", dpi=150)
    plt.close(fig)


def main():
    SONUC.mkdir(exist_ok=True)
    rng = np.random.default_rng(TOHUM)
    print(f"kontrol: G(B/2) = {kazanc_db(15.0, 30.0):.2f} dB (beklenen -3)")

    satirlar = []
    for B in ISIN_GENISLIKLERI:
        for s in SIGMALAR:
            satirlar += hata_olc(B, s, rng)
    df = pd.DataFrame(satirlar)
    df.to_csv(SONUC / "yon_hata.csv", index=False)
    print(f"gürültüsüz tepe RMSE beklenen: adım / sqrt(12) = {ADIM_DEG / np.sqrt(12):.2f} derece")
    print(df.pivot_table(index=["B_deg", "sigma_db"], columns="yontem", values="rmse_deg").round(2).to_string())
    grafikler(df, rng)
    print("yazıldı: results/yon_hata.csv, yon_hata.png, yon_ornek_tarama.png")


if __name__ == "__main__":
    main()
