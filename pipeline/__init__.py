"""Pipeline de analisis estadistico reproducible para Tomillo x Fusarium spp.

Paquete importable con un modulo por fase del flujo de trabajo:

- cargar_datos: carga y auditoria de calidad (fases 1-2)
- limpiar: normalizacion del dataset maestro (fase 3)
- eda: exploracion y descriptivos (fase 4)
- diseno: inferencia del diseno experimental (fase 5)
- supuestos: verificacion de supuestos de modelos (fase 6)
- modelos: analisis inferencial y seleccion automatica (fase 7)
- comparaciones: comparaciones multiples y letras CLD (fase 8)
- visualizar: figuras de resultados (fase 9)
- multivariado: PCA, clustering y categorias biologicas (fase 10)
- ranking: ranking de tecnicas de extraccion (fase 11)
- informe: generacion del informe final (fase 12)

El notebook orquestador ``dca/analisis_dca.ipynb`` importa este paquete.
"""

__version__ = "1.0.0"
