"""Adım 4: Hareketsiz hedef için RSSI'nın skaler Kalman filtresiyle süzülmesi.

Durum x: gerçek RSSI (dBm). Model x_k = x_(k-1) + w, w ~ N(0, Q). Ölçüm z_k = x_k + v, v ~ N(0, R), R = sigma^2.
Her adımda:
  tahmin:     P = P + Q
  kazanç:     K = P / (P + R)
  güncelleme: x = x + K (z - x),  P = (1 - K) P
Q = 0 ve P0 = R iken K_k = 1 / k olur ve filtre o ana kadarki ölçümlerin ortalamasına eşittir; main bunu
sayısal olarak da gösterir.

Not: Benzetimde gürültü her pakette bağımsız çekilir. Gerçekte hareketsiz bir hedefin gölgelemesi (engel kaybı)
paketten pakete aynı kalır; filtre yalnızca hızlı değişen bileşeni (sönümleme, ölçüm gürültüsü) azaltabilir.

Çıktılar (results/):
  kalman_ozet.csv       ham, hareketli ortalama ve Kalman için RSSI ve mesafe hatası
  kalman_yakinsama.csv  ölçüm sayısına göre RSSI RMSE
  kalman.png            tek bir koşu ve yakınsama grafiği
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ortak import SONUC, TOHUM, mesafe_kestir, rssi_uret

D_GERCEK = 100.0     # m
N_YOL = 2.7          # yol kaybı üsteli
SIGMA = 6.0          # dB, gölgeleme gürültüsü
OLCUM = 50           # paket sayısı
PENCERE = 10         # hareketli ortalama penceresi
Q_LISTESI = [0.0, 0.1]  # süreç gürültüsü (dB^2); 0 durağan hedef, 0.1 yavaş değişime izin verir
KOSU = 4000


def kalman(z, R, Q=0.0):
    """Skaler Kalman filtresi. z (koşu, ölçüm) boyutlu; her adımdaki kestirim ve P döner."""
    x = z[:, 0].copy()        # ilk ölçümle başla
    P = np.full(len(z), R)    # ilk belirsizlik ölçüm gürültüsü kadar
    xs, Ps = [x.copy()], [P.copy()]
    for k in range(1, z.shape[1]):
        P = P + Q
        K = P / (P + R)
        x = x + K * (z[:, k] - x)
        P = (1 - K) * P
        xs.append(x.copy())
        Ps.append(P.copy())
    return np.column_stack(xs), np.column_stack(Ps)


def hareketli_ortalama(z, pencere=PENCERE):
    """Son `pencere` ölçümün ortalaması; başta elde ne kadar ölçüm varsa onların ortalaması."""
    toplam = np.cumsum(z, axis=1)
    k = np.arange(1, z.shape[1] + 1)
    onceki = np.concatenate([np.zeros((len(z), pencere)), toplam[:, :-pencere]], axis=1)[:, : z.shape[1]]
    return (toplam - onceki) / np.minimum(k, pencere)


def main():
    SONUC.mkdir(exist_ok=True)
    rng = np.random.default_rng(TOHUM)
    gercek = rssi_uret(D_GERCEK, N_YOL, 0.0, rng)  # gürültüsüz RSSI, dBm
    z = rssi_uret(np.full((KOSU, OLCUM), D_GERCEK), N_YOL, SIGMA, rng)

    # Sağlama 1: gürültü yoksa filtre gerçek değerde kalır
    x0, _ = kalman(np.full((1, OLCUM), gercek), SIGMA**2)
    print(f"gürültüsüz Kalman en büyük sapma: {np.abs(x0 - gercek).max():.2e} dB")
    # Sağlama 2: Q = 0 iken Kalman = kümülatif ortalama
    xk, Pk = kalman(z, SIGMA**2, 0.0)
    kumulatif = np.cumsum(z, axis=1) / np.arange(1, OLCUM + 1)
    print(f"Q = 0 Kalman ile kümülatif ortalama arasındaki en büyük fark: {np.abs(xk - kumulatif).max():.2e} dB")
    print(f"Q = 0, {OLCUM} ölçüm sonunda P = {Pk[0, -1]:.4f}, beklenen R / N = {SIGMA**2 / OLCUM:.4f}")

    kestirimler = {"ham (tek paket)": z, f"hareketli ortalama ({PENCERE})": hareketli_ortalama(z)}
    for Q in Q_LISTESI:
        kestirimler[f"Kalman Q={Q:g}"] = kalman(z, SIGMA**2, Q)[0]

    satirlar, yakinsama = [], {"olcum_sayisi": np.arange(1, OLCUM + 1)}
    for ad, x in kestirimler.items():
        son = x[:, -1]
        bagil = np.abs(mesafe_kestir(son, N_YOL) - D_GERCEK) / D_GERCEK
        satirlar.append({
            "yontem": ad,
            "rssi_rmse_db": float(np.sqrt(np.mean((son - gercek) ** 2))),
            "mesafe_ort_bagil_hata": float(bagil.mean()),
            "mesafe_medyan_bagil_hata": float(np.median(bagil)),
            "mesafe_p90_bagil_hata": float(np.percentile(bagil, 90)),
        })
        yakinsama[ad] = np.sqrt(np.mean((x - gercek) ** 2, axis=0))
    ozet = pd.DataFrame(satirlar)
    ham = ozet.loc[0, "mesafe_ort_bagil_hata"]
    ozet["mesafe_hatasi_dusus"] = 1 - ozet["mesafe_ort_bagil_hata"] / ham
    ozet.to_csv(SONUC / "kalman_ozet.csv", index=False)
    pd.DataFrame(yakinsama).to_csv(SONUC / "kalman_yakinsama.csv", index=False)
    print(f"d = {D_GERCEK:g} m, n = {N_YOL}, sigma = {SIGMA:g} dB, {OLCUM} paket, {KOSU} koşu, "
          f"son paketten sonraki durum:")
    print(ozet.round(4).to_string(index=False))
    print(f"kuramsal: ham {SIGMA:.2f} dB, hareketli ortalama {SIGMA / np.sqrt(PENCERE):.2f} dB, "
          f"Kalman Q=0 {SIGMA / np.sqrt(OLCUM):.2f} dB")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    k = np.arange(1, OLCUM + 1)
    a1.plot(k, z[0], ".", color="0.6", label="raw RSSI")
    etiketler = {"ham (tek paket)": "raw (single packet)",
                 f"hareketli ortalama ({PENCERE})": f"moving average ({PENCERE})"}
    renk = dict(zip(kestirimler, ["0.5", "#ff7f0e", "#2ca02c", "#d62728"]))  # iki panelde aynı renk
    for ad, x in list(kestirimler.items())[1:]:
        a1.plot(k, x[0], "-", color=renk[ad], label=etiketler.get(ad, ad))
    a1.axhline(gercek, color="k", ls=":", label="true RSSI")
    a1.set_xlabel("packet")
    a1.set_ylabel("RSSI [dBm]")
    a1.set_title(f"One run, d = {D_GERCEK:g} m, σ = {SIGMA:g} dB")
    for ad in kestirimler:
        a2.plot(k, yakinsama[ad], color=renk[ad], label=etiketler.get(ad, ad))
    a2.plot(k, SIGMA / np.sqrt(k), "k:", label="σ / √k")
    a2.set_xlabel("packet")
    a2.set_ylabel("RSSI RMSE [dB]")
    a2.set_title(f"RMSE over {KOSU} runs")
    for ax in (a1, a2):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(SONUC / "kalman.png", dpi=150)
    plt.close(fig)
    print("yazıldı: results/kalman_ozet.csv, kalman_yakinsama.csv, kalman.png")


if __name__ == "__main__":
    main()
