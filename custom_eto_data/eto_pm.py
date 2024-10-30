import numpy as np

def calc_eto(df, mt, sr, ws, vps, hf=0, P=101.3, h=0.7):
    """
    Calcula la ETo para un Dataframe usando la fórmula de Penman-Monteith.
    Parámetros:
        df (dataframe): DataFrame que contiene los datos.
        mt (str)       : Nombre de la columna de Temperatura media (ºC).
        sr (str)       : Nombre de la columna de Radiación solar (W/m²).
        hf (str/int)   : Nombre de la columna de flujo de calor (si no existe, es 0).
        ws (str)       : Nombre de la columna de Velocidad del viento (m/s).
        vps (str)      : Nombre de la columna de Pendiente de presión de vapor (kPa/ºC).
        P (float)      : Presión atmosférica (kPa), dependiente de la altitud (por defecto 101.3).
        h (float)      : Humedad relativa, por defecto se fija en 70% (0.7).
    """

    # Constante psicrométrica (kPa/ºC)
    gamma = 0.000665 * P

    # Parámetros necesarios
    T = df[mt]
    Rn = df[sr] * 0.0864  # Conversión a MJ/m²/día con 24 horas de radiación
    G = df[hf] * 0.0864 if not isinstance(hf, int) else 0  # G = 0 si no se especifica flujo de calor
    u2 = df[ws]
    Delta = df[vps]

    # Presión de vapor de saturación (kPa)
    e_s = 0.6108 * np.exp((17.27 * T) / (T + 237.3))

    # Presión actual (usando la humedad relativa)
    e_a = e_s * h

    # Aplicar fórmula Penman-Monteith
    ETo = (0.408 * Delta * (Rn - G) + gamma * (900 / (T + 273)) * u2 * (e_s - e_a)) / (Delta + gamma * (1 + 0.34 * u2))

    # Añadir la columna ETo al DataFrame
    df['ETo (mm/día)'] = ETo

    return df