"""Adım 1: RSSI ile mesafe kestirimi ve hatanın n ve sigma'ya bağlılığı.

Çıktılar (results/):
  yol_kaybi_hata.csv    her mesafe ve (n, sigma) çifti için hata ölçütleri
  yol_kaybi_ozet.csv    (n, sigma) başına özet ve formülle karşılaştırma
  yol_kaybi_hata.png    bağıl hata - mesafe grafiği
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import erf

from ortak import F_HZ, SONUC, TOHUM, fspl_db, mesafe_kestir, rssi_d0, rssi_uret

MESAFELER = np.logspace(1, 3, 30)  # 10 ile 1000 m, logaritmik aralık
TEKRAR = 2000                      # her mesafe için gürültü tekrarı (PROJE: en az 500)
# Önce sigma sabitken n değişir, sonra n sabitken sigma değişir
SENARYOLAR = [(2.0, 6.0), (2.7, 6.0), (3.5, 6.0), (2.7, 3.0), (2.7, 8.0)]


def tarama(n, sigma_db, rng, mesafeler=MESAFELER, tekrar=TEKRAR):
    """Her mesafe için `tekrar` kez RSSI üretip mesafeyi geri hesaplar."""
    d = np.repeat(mesafeler[:, None], tekrar, axis=1)  # (mesafe, tekrar) boyutlu
    d_kestirim = mesafe_kestir(rssi_uret(d, n, sigma_db, rng), n)
    bagil = np.abs(d_kestirim - d) / d
    return pd.DataFrame({
        "d_m": mesafeler,
        "n": n,
        "sigma_db": sigma_db,
        "ort_bagil_hata": bagil.mean(axis=1),
        "p90_bagil_hata": np.percentile(bagil, 90, axis=1),
        "std_ln_oran": np.log(d_kestirim / d).std(axis=1),
        "medyan_oran": np.median(d_kestirim / d, axis=1),
        "ortalama_oran": (d_kestirim / d).mean(axis=1),
    })


def grafik(df):
    fig, eksenler = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    gruplar = [("n changes, sigma = 6 dB", df["sigma_db"] == 6.0),
               ("sigma changes, n = 2.7", df["n"] == 2.7)]
    for ax, (baslik, secim) in zip(eksenler, gruplar):
        for (n, s), k in df[secim].groupby(["n", "sigma_db"]):
            cizgi, = ax.plot(k["d_m"], 100 * k["ort_bagil_hata"], "-", label=f"n={n}, σ={s:g} dB mean")
            ax.plot(k["d_m"], 100 * k["p90_bagil_hata"], "--", color=cizgi.get_color(),
                    label=f"n={n}, σ={s:g} dB 90%")
        ax.set_xscale("log")
        ax.set_xlabel("true distance [m]")
        ax.set_title(baslik)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    eksenler[0].set_ylabel("relative distance error [%]")
    fig.tight_layout()
    fig.savefig(SONUC / "yol_kaybi_hata.png", dpi=150)
    plt.close(fig)


def main():
    SONUC.mkdir(exist_ok=True)
    rng = np.random.default_rng(TOHUM)

    print(f"FSPL 868 MHz, 100 m: {fspl_db(100, 868e6):.2f} dB  (beklenen yaklaşık 71)")
    print(f"FSPL 433 MHz, 100 m: {fspl_db(100, 433e6):.2f} dB  (beklenen yaklaşık 65)")
    print(f"RSSI(d0 = 1 m) = {rssi_d0():.2f} dBm  (f = {F_HZ / 1e6:.0f} MHz, P_tx = 14 dBm)")

    gurultusuz = tarama(2.7, 0.0, rng, tekrar=1)
    print(f"gürültüsüz en büyük bağıl hata: {gurultusuz['ort_bagil_hata'].max():.2e}")

    df = pd.concat([tarama(n, s, rng) for n, s in SENARYOLAR], ignore_index=True)
    df.to_csv(SONUC / "yol_kaybi_hata.csv", index=False)

    # Hata mesafeden bağımsız olmalı: d_kestirim / d = 10^(X / (10 n))
    ozet = df.groupby(["n", "sigma_db"]).agg(
        ort_bagil_hata=("ort_bagil_hata", "mean"),
        p90_bagil_hata=("p90_bagil_hata", "mean"),
        std_ln_oran=("std_ln_oran", "mean"),
        ort_hata_min=("ort_bagil_hata", "min"),
        ort_hata_max=("ort_bagil_hata", "max"),
        medyan_oran=("medyan_oran", "mean"),
        ortalama_oran=("ortalama_oran", "mean"),
    ).reset_index()
    ozet["formul_ln10_10n_sigma"] = np.log(10) / (10 * ozet["n"]) * ozet["sigma_db"]
    # d_kestirim / d = exp(s Z), Z ~ N(0, 1). Bu log-normal dağılımın kapalı formları:
    #   std(ln oran) = s,  medyan oran = 1,  ortalama oran = exp(s^2 / 2),
    #   E|oran - 1| = exp(s^2 / 2) * erf(s / sqrt(2))
    s = ozet["formul_ln10_10n_sigma"]
    ozet["beklenen_ort_bagil_hata"] = np.exp(s**2 / 2) * erf(s / np.sqrt(2))
    ozet["beklenen_ortalama_oran"] = np.exp(s**2 / 2)
    ozet.to_csv(SONUC / "yol_kaybi_ozet.csv", index=False)
    print(ozet.round(4).to_string(index=False))

    grafik(df)
    print("yazıldı: results/yol_kaybi_hata.csv, yol_kaybi_ozet.csv, yol_kaybi_hata.png")


if __name__ == "__main__":
    main()
