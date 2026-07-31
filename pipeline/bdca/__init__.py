"""Pipeline de análisis de rendimiento de BDCA (Jenkyn mildiu).

Paquete importable con un módulo por fase del flujo de trabajo para el diseño
RCBD (bloque aleatorizado) de una sola respuesta (rendimiento):

- cargar: carga y auditoría de calidad del CSV BDCA (fases 1-2)
- eda: exploración descriptiva por tratamiento y distribución (fase 4)
- supuestos: verificación de supuestos del modelo de bloques (fase 6)
- modelos: análisis inferencial de bloques (ANOVA clásico + LMM + ICC) (fase 7)
- comparaciones: comparación múltiple post-hoc (Tukey HSD, pares vs referencia R) (fase 8)
- informe: generación del informe final con resumen del diseño RCBD (fase 12)

El notebook orquestador ``bdca/analisis_bdca.ipynb`` importa este paquete.
"""

__version__ = "1.0.0"

from pipeline.bdca.cargar import cargar
from pipeline.bdca.eda import resumen_descriptivo, figuras_eda
from pipeline.bdca.supuestos import analisis_supuestos
from pipeline.bdca.modelos import anova_bloques, lmm_bloques
from pipeline.bdca.comparaciones import posthoc_tukey
from pipeline.bdca.informe import informe_completo
