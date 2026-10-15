"""Adım 3: İki antenli interferometre ile faz farkından geliş açısı.

Faz farkı:  delta_phi = 2 pi d sin(theta) / lambda   (radyan; d anten aralığı m, theta derece)
Kestirim:   theta = arcsin(lambda delta_phi / (2 pi d))
Faz [-pi, pi) aralığında ölçülür. d > lambda / 2 olursa birden çok theta aynı ölçülen faza düşer.

Çıktılar (results/):
  faz_belirsizlik.png   d = lambda/2 ve d = 2 lambda için gerçek açı - kestirilen açı (gürültüsüz)
  faz_belirsizlik.csv   anten aralığına göre aday açı sayısı ve belirsizliksiz bölge
  faz_hata.csv          açı ve faz gürültüsüne göre hata, formülle karşılaştırma
  faz_hata.png          hata - geliş açısı
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ortak import F_HZ, SONUC, TOHUM, dalga_boyu

LAMBDA = dalga_boyu(F_HZ)              # 868 MHz için 0,3456 m
SIGMA_PHI_DEG = [2.0, 5.0, 10.0]       # faz ölçüm gürültüsü (derece)
TETALAR = np.array([0, 10, 20, 30, 40, 50, 60, 70, 75, 80, 85, 88])  # derece
TEKRAR = 5000


def sar(faz_rad):
    """Fazı [-pi, pi) aralığına sarar (radyan)."""
    return (faz_rad + np.pi) % (2 * np.pi) - np.pi


def faz_olc(theta_deg, d_m, sigma_phi_rad, rng, lam=LAMBDA):
    """Ölçülen faz farkı (radyan): gerçek faz + gürültü, sonra sarma."""
    gercek = 2 * np.pi * d_m * np.sin(np.radians(theta_deg)) / lam
    return sar(gercek + rng.normal(0.0, sigma_phi_rad, np.shape(theta_deg)))


def aci_kestir(dphi_rad, d_m, lam=LAMBDA):
    """Ana (k = 0) çözüm, derece. Gürültü arcsin argümanını 1'in üstüne itebilir, kırpılır."""
    return np.degrees(np.arcsin(np.clip(lam * dphi_rad / (2 * np.pi * d_m), -1.0, 1.0)))


def adaylar(dphi_rad, d_m, lam=LAMBDA):
    """Ölçülen faza uyan bütün geliş açıları: delta_phi + 2 pi k, |sin theta| <= 1 olan her k."""
    kmax = int(np.ceil(d_m / lam)) + 1
    s = lam * (dphi_rad + 2 * np.pi * np.arange(-kmax, kmax + 1)) / (2 * np.pi * d_m)
    return np.degrees(np.arcsin(s[np.abs(s) <= 1.0]))


def belirsizlik(rng):
    tum = np.linspace(-89.5, 89.5, 359)
    fig, eksenler = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, oran in zip(eksenler, [0.5, 2.0]):
        d = oran * LAMBDA
        dphi = faz_olc(tum, d, 0.0, rng)
        for t, f in zip(tum[::6], dphi[::6]):
            a = adaylar(f, d)
            ax.plot(np.full(len(a), t), a, ".", color="0.7", ms=3)
        ax.plot(tum, aci_kestir(dphi, d), "-", color="#1f77b4", lw=2, label="principal solution (k = 0)")
        ax.plot([-90, 90], [-90, 90], "k:", lw=1, label="ideal")
        ax.plot([], [], ".", color="0.7", label="all angles with the same phase")
        if oran > 0.5:
            sinir = np.degrees(np.arcsin(1 / (2 * oran)))
            ax.axvspan(-sinir, sinir, color="#2ca02c", alpha=0.12, label=f"unambiguous |θ| < {sinir:.1f}°")
        ax.set_title(f"d = {oran:g} λ = {d * 100:.1f} cm")
        ax.set_xlabel("true angle of arrival θ [deg]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
    eksenler[0].set_ylabel("estimated angle [deg]")
    fig.tight_layout()
    fig.savefig(SONUC / "faz_belirsizlik.png", dpi=150)
    plt.close(fig)

    # Anten aralığına göre aday sayısı: yaklaşık 2 d / lambda, belirsizliksiz bölge 2 arcsin(lambda / 2d)
    satirlar = []
    for oran in [0.25, 0.5, 1.0, 2.0, 4.0]:
        d = oran * LAMBDA
        sayi = [len(adaylar(f, d)) for f in faz_olc(tum, d, 0.0, rng)]
        satirlar.append({
            "d_bolu_lambda": oran, "d_cm": 100 * d,
            "ort_aday_sayisi": float(np.mean(sayi)), "en_cok_aday": int(np.max(sayi)),
            "belirsizliksiz_bolge_deg": float(2 * np.degrees(np.arcsin(min(1.0, 1 / (2 * oran))))),
        })
    return pd.DataFrame(satirlar)


def hata(rng):
    satirlar = []
    for oran in [0.5, 2.0]:
        d = oran * LAMBDA
        for s_deg in SIGMA_PHI_DEG:
            s = np.radians(s_deg)
            for t in TETALAR:
                teta = np.full(TEKRAR, float(t))
                dphi = faz_olc(teta, d, s, rng)
                if oran <= 0.5:
                    kestirim = aci_kestir(dphi, d)
                else:
                    # d = 2 lambda: belirsizlik çözülmüş kabul edilir (doğru k biliniyor).
                    # Yalnızca gürültünün etkisini görmek için; gerçekte k bilinmez.
                    gercek_faz = 2 * np.pi * d * np.sin(np.radians(teta)) / LAMBDA
                    acik = gercek_faz + sar(dphi - gercek_faz)
                    kestirim = np.degrees(np.arcsin(np.clip(LAMBDA * acik / (2 * np.pi * d), -1, 1)))
                e = np.abs(kestirim - teta)
                satirlar.append({
                    "d_bolu_lambda": oran, "sigma_phi_deg": s_deg, "theta_deg": t,
                    "rmse_deg": float(np.sqrt(np.mean(e**2))),
                    "medyan_deg": float(np.median(e)),
                    "p90_deg": float(np.percentile(e, 90)),
                    # Formül: sigma_theta = sigma_phi lambda / (2 pi d cos(theta)), küçük gürültü yaklaşımı
                    "formul_deg": float(np.degrees(s * LAMBDA / (2 * np.pi * d * np.cos(np.radians(t))))),
                    "ters_isaret_orani": float(np.mean(np.sign(kestirim) != np.sign(teta))) if t > 0 else np.nan,
                })
    return pd.DataFrame(satirlar)


def hata_grafigi(df):
    fig, eksenler = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, oran in zip(eksenler, [0.5, 2.0]):
        for s_deg in SIGMA_PHI_DEG:
            k = df[(df["d_bolu_lambda"] == oran) & (df["sigma_phi_deg"] == s_deg)]
            cizgi, = ax.plot(k["theta_deg"], k["rmse_deg"], "o-", label=f"simulated RMSE, σφ = {s_deg:g}°")
            ax.plot(k["theta_deg"], k["formul_deg"], "--", color=cizgi.get_color(),
                    label=f"σφ λ / (2π d cos θ), σφ = {s_deg:g}°")
        baslik = f"d = {oran:g} λ" + ("" if oran <= 0.5 else " (correct k assumed)")
        ax.set_title(baslik)
        ax.set_yscale("log")
        ax.set_xlabel("true angle of arrival θ [deg]")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
    eksenler[0].set_ylabel("angle error [deg]")
    fig.tight_layout()
    fig.savefig(SONUC / "faz_hata.png", dpi=150)
    plt.close(fig)


def main():
    SONUC.mkdir(exist_ok=True)
    rng = np.random.default_rng(TOHUM)
    print(f"lambda = {LAMBDA:.4f} m (868 MHz), lambda/2 = {100 * LAMBDA / 2:.2f} cm")

    tum = np.linspace(-89, 89, 179)
    for oran in [0.5, 2.0]:
        d = oran * LAMBDA
        e = np.abs(aci_kestir(faz_olc(tum, d, 0.0, rng), d) - tum)
        print(f"gürültüsüz, d = {oran:g} lambda: en büyük hata {e.max():.2e} derece, "
              f"hatası 1 dereceden büyük açı oranı {np.mean(e > 1):.1%}")
    print("örnek: theta = 40 derece, d = 2 lambda için aynı faza düşen açılar:",
          np.round(adaylar(faz_olc(np.array(40.0), 2 * LAMBDA, 0.0, rng), 2 * LAMBDA), 2))

    tablo = belirsizlik(rng)
    tablo.to_csv(SONUC / "faz_belirsizlik.csv", index=False)
    print(tablo.round(2).to_string(index=False))

    df = hata(rng)
    df.to_csv(SONUC / "faz_hata.csv", index=False)
    print(df[df["sigma_phi_deg"] == 5.0][["d_bolu_lambda", "theta_deg", "rmse_deg", "medyan_deg",
                                          "formul_deg", "ters_isaret_orani"]].round(3).to_string(index=False))
    hata_grafigi(df)
    print("yazıldı: results/faz_belirsizlik.png, faz_belirsizlik.csv, faz_hata.csv, faz_hata.png")


if __name__ == "__main__":
    main()
