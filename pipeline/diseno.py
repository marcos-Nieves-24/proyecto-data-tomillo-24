"""Fase 5: inferencia del diseño experimental.

Documenta el diseño subyacente (DCA factorial método × aislado, 3 réplicas
biológicas), la unidad experimental y las limitaciones (control compartido).
"""

from __future__ import annotations

import pandas as pd

from pipeline.config import (
    CONCENTRACION_UNICA_MG_ML,
    VARIABLES_RESPUESTA,
    VARIABLE_LABEL,
    guardar_tabla,
)


def inferir_diseno(df_bio: pd.DataFrame) -> dict:
    """Infiera y documenta el diseño experimental del bioensayo.

    Devuelve un dict con el detalle tabular, atributos clave y un texto
    explicativo que justifica el diseño y sus límites.
    """
    n_aislados = df_bio["aislamiento"].nunique()
    n_metodos = df_bio["metodo_extraccion"].nunique()
    n_replicas = df_bio["replica"].nunique()
    n_tratamientos = int(n_aislados * n_metodos)
    conteo = df_bio.groupby(["metodo_extraccion", "aislamiento"]).size()
    balanceado = bool((conteo == 3).all())
    replicas_por_celda = int(conteo.iloc[0]) if len(conteo) else 0

    detalle = pd.DataFrame({
        "atributo": [
            "Tipo de diseño",
            "Factores",
            "Niveles de método",
            "Niveles de aislado",
            "Número de tratamientos (método × aislado)",
            "Réplicas biológicas por celda",
            "Unidad experimental",
            "Concentración ensayada (mg/mL)",
            "Variables de respuesta",
            "Diseño balanceado",
        ],
        "valor": [
            "DCA factorial método × aislado",
            "método de extracción (fijo); aislado (fijo en ANOVA, aleatorio en LMM)",
            str(n_metodos),
            str(n_aislados),
            str(n_tratamientos),
            str(replicas_por_celda),
            "Caja Petri",
            str(CONCENTRACION_UNICA_MG_ML),
            "; ".join(VARIABLE_LABEL[v] for v in VARIABLES_RESPUESTA),
            "Sí" if balanceado else "No",
        ],
    })
    guardar_tabla(detalle, "diseno_experimental", index=False)

    texto = (
        "El diseño se infiere como un DCA factorial con dos factores fijos "
        "(técnica de extracción × aislado de Fusarium) y tres réplicas "
        "biológicas por combinación; la unidad experimental es la caja Petri. "
        "La concentración (5 mg/mL) es constante y NO es factor experimental. "
        "Caveat: cada porcentaje de inhibición se calculó contra un único "
        "control C4 compartido por las tres réplicas del aislado; por lo tanto, "
        "las réplicas de %INH no son totalmente independientes "
        "(pseudorreplicación del control). El análisis del crecimiento crudo "
        "(mm) no presenta este problema."
    )

    print("=" * 72)
    print("FASE 5 - Diseño experimental inferido")
    print("=" * 72)
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(detalle.to_string(index=False))
    print(f"\n{texto}")

    return {
        "n_tratamientos": n_tratamientos,
        "n_aislados": int(n_aislados),
        "n_metodos": int(n_metodos),
        "n_replicas": int(n_replicas),
        "replicas_por_celda": replicas_por_celda,
        "balanceado": balanceado,
        "detalle": detalle,
        "texto": texto,
    }
