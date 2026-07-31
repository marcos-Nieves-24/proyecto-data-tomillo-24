"""
08_exportar.py — Exportación de datasets finales a resultados/tablas/.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import json
from config import DIR_DATOS, DIR_TABLAS

ARCHIVOS = [
    ("rendimiento_extraccion.csv", "rendimiento_extraccion.csv"),
    ("crecimiento_micelial_con_inhibicion.csv", "crecimiento_micelial.csv"),
    ("crecimiento_micelial_raw.csv", "crecimiento_micelial_raw.csv"),
    ("crecimiento_micelial_estandarizado.csv", "crecimiento_micelial_estandarizado.csv"),
    ("conidias_con_inhibicion.csv", "conidias.csv"),
    ("conidias_raw.csv", "conidias_raw.csv"),
    ("conidias_estandarizado.csv", "conidias_estandarizado.csv"),
]

print("📄 Exportando datasets finales...")
exportados = []
for origen, destino in ARCHIVOS:
    ruta_origen = DIR_DATOS / origen
    ruta_destino = DIR_TABLAS / destino
    if not ruta_origen.exists():
        print(f"  ⚠ No encontrado: {origen}")
        continue
    df = pd.read_csv(ruta_origen)
    df.to_csv(ruta_destino, index=False)
    exportados.append({"archivo": destino, "filas": len(df), "columnas": len(df.columns)})
    print(f"  ✅ {destino} ({len(df)} filas × {len(df.columns)} cols)")

meta = {
    "proyecto": "Tomillo × Fusarium spp.",
    "fecha": "2026-07-28",
    "nota": "Controles desnormalizados por fila. Ver docs/05_DATA_QUALITY_ISSUES.md",
    "tablas": exportados,
}
with open(DIR_TABLAS / "metadatos.json", "w", encoding="utf-8") as f:
    json.dump(meta, f, ensure_ascii=False, indent=2)

print(f"\n✅ Exportación completa. {len(exportados)} tablas exportadas.")
