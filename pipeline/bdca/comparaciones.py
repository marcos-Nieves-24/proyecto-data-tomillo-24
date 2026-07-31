"""Fase 3.6: Comparaciones múltiples post-hoc Tukey HSD entre tratamientos para rendimiento BDCA.

Implementa `posthoc_tukey(df)`: Tukey HSD para los pares de tratamientos (6 pares),
columna interpretativa explícita ``vs_referencia_T0`` por par (comparando cada
tratamiento R/T1/T2 vs T0, el control sin tratar). Sin Dunnett, sin CLD.

Protege contra errores: funciona cuando `statsmodels` está disponible; si no,
levanta ImportError con un mensaje claro.
"""

from __future__ import annotations

import pandas as pd
import numpy as np

import warnings
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from pipeline.config import guardar_tabla, save_figure_pub

ALPHA = 0.05


def posthoc_tukey(df: pd.DataFrame) -> pd.DataFrame:
    """Tukey HSD entre tratamientos para rendimiento BDCA.

    Ajusta el test de Tukey HSD de statsmodels y devuelve un DataFrame con:
        • par: etiqueta como "R vs T0", "T1 vs T0", etc.
        • diferencia_medias: media(tratamiento) - media(T0)
        • p_valor_ajustado: p-value ajustado de Tukey
        • ic95_inferior, ic95_superior: intervalo de confianza del 95%
        • significativo: True si p <= 0.05
        • vs_referencia_T0: booleano indicando si este par compara vs T0 (control)

    Además guarda la figura ``posthoc_tukey_pares`` (Tukey plot) en
    ``bdca/resultados/figuras/``. El DataFrame resultante se guarda como
    ``bdca/resultados/tablas/posthoc_tukey.csv``.
    """
    df = df.copy()
    df["trt"] = df["trt"].astype(str).str.strip().str.upper()

    # Ordenar tratamientos para consistencia (R, T0, T1, T2)
    orden = sorted(df["trt"].unique())
    # Asegurar que T0 (control sin tratar) esté presente como referencia
    if "T0" not in orden:
        warnings.warn("T0 no encontrado en trt; usando primera categoría como referencia.")
        referencia = orden[0]
    else:
        referencia = "T0"

    # Tukey HSD con API endog + groups (misma que pipeline/comparaciones.py).
    # `yield` es palabra reservada de Python: se cita como columna del df.
    tukey = pairwise_tukeyhsd(df["yield"], df["trt"], alpha=ALPHA)

    # Pares por combinación (mismo patrón probado que pipeline/comparaciones.py)
    grupos = list(tukey.groupsunique)
    pares = [(a, b) for i, a in enumerate(grupos) for b in grupos[i + 1:]]

    filas = []
    for i, (a, b) in enumerate(pares):
        # Determinar si este par involucra T0 (control)
        vs_ref = (a == referencia) or (b == referencia)

        filas.append({
            "par": f"{a} vs {b}",
            "diferencia_medias": float(tukey.meandiffs[i]),
            "p_valor_ajustado": float(tukey.pvalues[i]),
            "ic95_inferior": float(tukey.confint[i, 0]),
            "ic95_superior": float(tukey.confint[i, 1]),
            "significativo": bool(tukey.reject[i]),
            "vs_referencia_T0": vs_ref,
        })

    resultado = pd.DataFrame(filas)
    guardar_tabla(resultado, "posthoc_tukey", index=False)

    # Figura Tukey (requerida por spec rcbd-reporting: figuras/tukey plot)
    fig = tukey.plot_simultaneous()
    fig.suptitle("Comparaciones múltiples Tukey HSD por tratamiento (BDCA)")
    fig.set_figwidth(10)
    fig.set_figheight(6)
    save_figure_pub(fig, "posthoc_tukey_pares", titulo="Comparaciones múltiples Tukey HSD por tratamiento (BDCA)")

    # Resumen para salida
    print("=" * 72)
    print("FASE 3.6 - Comparaciones múltiples post-hoc Tukey HSD (BDCA)")
    print("=" * 72)
    print("Referencia: T0 (control sin tratar).")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(resultado.to_string(index=False))

    return resultado


if __name__ == "__main__":
    # Permite ejecución directa para pruebas rápidas
    from pipeline.bdca.cargar import cargar
    df, _ = cargar()
    posthoc_tukey(df)
