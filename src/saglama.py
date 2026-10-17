"""PROJE.md bölüm 7'deki sağlama kontrolleri, tek komutla.

Çalıştırma: python src/saglama.py
Her kontrol assert ile sınanır; biri tutmazsa betik hata verip durur.
"""
import numpy as np

import faz_farki as faz
import kalman as kal
import trilaterasyon as tri
import yol_kaybi as yk
import yonlu_anten as ya
from ortak import TOHUM, dalga_boyu, fspl_db, mesafe_kestir, rssi_uret

sayac = 0


def kontrol(ad, kosul, ayrinti):
    global sayac
    assert kosul, f"TUTMADI: {ad} ({ayrinti})"
    sayac += 1
    print(f"[tamam] {ad}: {ayrinti}")


rng = np.random.default_rng(TOHUM)

# Temel bilgi (PROJE.md bölüm 5)
for f, beklenen in ((433e6, 0.69), (868e6, 0.345), (121.5e6, 2.47)):
    kontrol(f"dalga boyu {f / 1e6:g} MHz", abs(dalga_boyu(f) - beklenen) < 0.01, f"{dalga_boyu(f):.4f} m")
kontrol("FSPL 868 MHz 100 m", abs(fspl_db(100, 868e6) - 71) < 0.5, f"{fspl_db(100, 868e6):.2f} dB")
kontrol("FSPL 433 MHz 100 m", abs(fspl_db(100, 433e6) - 65) < 0.5, f"{fspl_db(100, 433e6):.2f} dB")

# Adım 1: gürültüsüz mesafe tam, gürültülü hata formülle uyumlu
d = np.logspace(1, 3, 50)
e = np.abs(mesafe_kestir(rssi_uret(d, 2.7, 0.0, rng), 2.7) - d) / d
kontrol("gürültüsüz mesafe", e.max() < 1e-9, f"en büyük bağıl hata {e.max():.1e}")
t = yk.tarama(2.7, 6.0, rng)
s = np.log(10) / (10 * 2.7) * 6.0
kontrol("std(ln d̂/d) = ln10 / (10 n) σ", abs(t["std_ln_oran"].mean() / s - 1) < 0.02,
        f"benzetim {t['std_ln_oran'].mean():.4f}, formül {s:.4f}")
kontrol("hata mesafeden bağımsız", np.ptp(t["std_ln_oran"]) < 0.06,
        f"30 mesafede std(ln) {t['std_ln_oran'].min():.3f} ile {t['std_ln_oran'].max():.3f} arası")

# Adım 2: kazanç tanımı ve gürültüsüz yön
kontrol("G(B/2) = -3 dB", abs(ya.kazanc_db(15.0, 30.0) + 3) < 1e-12, f"{ya.kazanc_db(15.0, 30.0):.3f} dB")
gercek = rng.uniform(0, 360, 2000)
for B in ya.ISIN_GENISLIKLERI:
    rssi = ya.tarama_uret(gercek, B, 0.0, rng)
    e_parabol = np.abs(ya.aci_farki(ya.kestir_parabol(rssi), gercek)).max()
    e_gauss = np.abs(ya.aci_farki(ya.kestir_gauss(rssi, B), gercek)).max()
    e_tepe = np.abs(ya.aci_farki(ya.kestir_tepe(rssi), gercek)).max()
    kontrol(f"gürültüsüz yön, B = {B:g}", max(e_parabol, e_gauss) < 1e-6 and e_tepe <= ya.ADIM_DEG / 2 + 1e-9,
            f"parabol {e_parabol:.1e}, gauss {e_gauss:.1e}, tepe en çok {e_tepe:.2f} (≤ adım/2)")

# Adım 3: faz belirsizliği ve gürültü formülü
teta = np.linspace(-89, 89, 179)
for oran in (0.25, 0.5):
    dd = oran * faz.LAMBDA
    e = np.abs(faz.aci_kestir(faz.faz_olc(teta, dd, 0.0, rng), dd) - teta).max()
    kontrol(f"d = {oran:g} λ belirsiz değil", e < 1e-9, f"gürültüsüz en büyük hata {e:.1e} derece")
dd = 2 * faz.LAMBDA
e = np.abs(faz.aci_kestir(faz.faz_olc(teta, dd, 0.0, rng), dd) - teta)
sinir = np.degrees(np.arcsin(faz.LAMBDA / (2 * dd)))
kontrol("d = 2 λ belirsiz", (e[np.abs(teta) < sinir - 1] < 1e-9).all() and (e[np.abs(teta) > sinir + 1] > 1).all(),
        f"|θ| < {sinir:.1f}° doğru, dışında hepsi yanlış ({np.mean(e > 1):.0%} yanlış)")
aday = faz.adaylar(faz.faz_olc(np.array(40.0), dd, 0.0, rng), dd)
kontrol("d = 2 λ aday sayısı ≈ 2d/λ", len(aday) == 4 and np.isclose(aday, 40.0).any(), f"{np.round(aday, 2)}")
sp = np.radians(5.0)
for t_deg in (0, 30, 60):
    e = faz.aci_kestir(faz.faz_olc(np.full(20000, float(t_deg)), faz.LAMBDA / 2, sp, rng), faz.LAMBDA / 2) - t_deg
    beklenen = np.degrees(sp / (np.pi * np.cos(np.radians(t_deg))))
    kontrol(f"faz gürültüsü θ = {t_deg}°", abs(e.std() / beklenen - 1) < 0.05,
            f"benzetim {e.std():.3f}°, formül {beklenen:.3f}°")

# Adım 4: Kalman
z = rssi_uret(np.full((500, 50), 100.0), 2.7, 6.0, rng)
x, P = kal.kalman(z, 36.0, 0.0)
fark = np.abs(x - np.cumsum(z, axis=1) / np.arange(1, 51)).max()
kontrol("Kalman Q=0 = kümülatif ortalama", fark < 1e-9, f"en büyük fark {fark:.1e} dB")
kontrol("Kalman P_N = R / N", np.isclose(P[0, -1], 36.0 / 50), f"P = {P[0, -1]:.4f}")
gercek_rssi = rssi_uret(100.0, 2.7, 0.0, rng)
x0, _ = kal.kalman(np.full((1, 50), gercek_rssi), 36.0)
kontrol("gürültüsüz Kalman", np.abs(x0 - gercek_rssi).max() < 1e-12, "gerçek değerde kalıyor")

# Adım 5: trilaterasyon
for ad in ("ucgen", "neredeyse_dogrusal"):
    a = tri.GEOMETRILER[ad]
    p = tri.konum_bul(a, np.linalg.norm(a - tri.HEDEF, axis=1), rng)
    kontrol(f"gürültüsüz konum, {ad}", np.linalg.norm(p - tri.HEDEF) < 1e-6,
            f"hata {np.linalg.norm(p - tri.HEDEF):.1e} m")
a = tri.GEOMETRILER["dogrusal"]
ayna = np.array([tri.HEDEF[0], -tri.HEDEF[1]])
bulunan = [tri.konum_bul(a, np.linalg.norm(a - tri.HEDEF, axis=1), rng) for _ in range(40)]
gercege = sum(np.linalg.norm(p - tri.HEDEF) < 1e-6 for p in bulunan)
aynaya = sum(np.linalg.norm(p - ayna) < 1e-6 for p in bulunan)
kontrol("doğrusal geometride ayna belirsizliği", gercege + aynaya == 40 and gercege > 0 and aynaya > 0,
        f"40 çözümün {gercege} tanesi gerçek, {aynaya} tanesi ayna noktası")

# Tekrar üretilebilirlik: aynı tohum aynı sayıyı verir
r1 = yk.tarama(2.7, 6.0, np.random.default_rng(TOHUM))
r2 = yk.tarama(2.7, 6.0, np.random.default_rng(TOHUM))
kontrol("tohum sabit", r1.equals(r2), "aynı tohumla iki koşu birebir aynı")

print(f"\n{sayac} kontrolün hepsi tuttu.")
