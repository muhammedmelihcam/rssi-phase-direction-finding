"""Ortak sabitler, yollar ve yol kaybı fonksiyonları.

Birimler: frekans Hz, uzunluk m, güç dBm, kayıp ve gürültü dB, açı derece (aksi yazmadıkça).
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # ekran açmadan PNG yazmak için
import numpy as np

C = 3e8          # ışık hızı, m/s
TOHUM = 0        # bütün benzetimler bu tohumla başlar, sonuçlar tekrar üretilebilir
F_HZ = 868e6     # varsayılan frekans (LoRa, Avrupa bandı)
P_TX_DBM = 14.0  # verici gücü, LoRa için tipik üst değer
D0_M = 1.0       # log-distance modelinin referans mesafesi

KOK = Path(__file__).resolve().parent.parent
SONUC = KOK / "results"
VERI = KOK / "data"


def dalga_boyu(f_hz):
    """Dalga boyu (m). lambda = c / f."""
    return C / f_hz


def fspl_db(d_m, f_hz):
    """Serbest uzay yol kaybı (dB). Formül km ve MHz ister, burada çevrilir."""
    return 20 * np.log10(np.asarray(d_m) / 1e3) + 20 * np.log10(f_hz / 1e6) + 32.44


def rssi_d0(f_hz=F_HZ, p_tx_dbm=P_TX_DBM, d0_m=D0_M):
    """Referans mesafedeki RSSI (dBm): verici gücü eksi d0'daki serbest uzay kaybı."""
    return p_tx_dbm - fspl_db(d0_m, f_hz)


def rssi_uret(d_m, n, sigma_db, rng, rssi0_dbm=None, d0_m=D0_M):
    """Log-distance modeliyle gürültülü RSSI (dBm).

    RSSI(d) = RSSI(d0) - 10 n log10(d / d0) + X,  X ~ N(0, sigma_db) gölgeleme gürültüsü.
    """
    rssi0_dbm = rssi_d0() if rssi0_dbm is None else rssi0_dbm
    d_m = np.asarray(d_m, dtype=float)
    return rssi0_dbm - 10 * n * np.log10(d_m / d0_m) + rng.normal(0.0, sigma_db, d_m.shape)


def mesafe_kestir(rssi_dbm, n, rssi0_dbm=None, d0_m=D0_M):
    """RSSI'dan mesafe (m), log-distance modelinin tersi."""
    rssi0_dbm = rssi_d0() if rssi0_dbm is None else rssi0_dbm
    return d0_m * 10 ** ((rssi0_dbm - np.asarray(rssi_dbm)) / (10 * n))


def aci_farki(a_deg, b_deg):
    """İki açı arasındaki işaretli fark, [-180, 180) aralığında (derece)."""
    return (np.asarray(a_deg) - np.asarray(b_deg) + 180.0) % 360.0 - 180.0
