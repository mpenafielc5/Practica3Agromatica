from Bio import SeqIO
from Bio.Seq import Seq
import csv

INPUT_FASTA = "ls_orchid.fasta.txt"
MOTIVOS = ["ATG", "GCGG", "TATA", "AAGG"]  # modificar según necesidad

# Estrategia para manejar codón parcial: "TRIM" o "PAD"
STRATEGY = "TRIM"

# Límite de posiciones a mostrar en el .txt por motivo (evita líneas enormes)
MAX_POS_SHOW = 10

def ajustar_longitud_para_traduccion(seq_str, strategy="TRIM"):
    """Devuelve secuencia ajustada (string). TRIM recorta, PAD añade 'N's."""
    mod = len(seq_str) % 3
    if mod == 0:
        return seq_str
    if strategy == "TRIM":
        return seq_str[: len(seq_str) - mod]
    elif strategy == "PAD":
        return seq_str + ("N" * (3 - mod))
    else:
        raise ValueError("STRATEGY debe ser 'TRIM' o 'PAD'")

def encontrar_posiciones(seq, motivo):
    """Retorna lista de posiciones 1-based donde aparece motivo (permite solapamiento)."""
    posiciones = []
    start = 0
    while True:
        pos = seq.find(motivo, start)
        if pos == -1:
            break
        posiciones.append(pos + 1)  # 1-based
        start = pos + 1
    return posiciones

resultados = []

for record in SeqIO.parse(INPUT_FASTA, "fasta"):
    seq_id = record.id
    seq_str = str(record.seq).upper()
    longitud = len(seq_str)

    # Ajustar secuencia para evitar codones parciales y la advertencia
    seq_para_traducir = ajustar_longitud_para_traduccion(seq_str, STRATEGY)

    # Traducción segura (si tiene al menos 3 nt después del ajuste)
    if len(seq_para_traducir) >= 3:
        try:
            traduccion = str(Seq(seq_para_traducir).translate(to_stop=False))
        except Exception as e:
            traduccion = f"ERROR_TRAD: {e}"
    else:
        traduccion = ""

    # Buscar motivos y registrar posiciones completas
    detalle_motivos_full = {}
    for motivo in MOTIVOS:
        pos_list = encontrar_posiciones(seq_str, motivo)  # buscamos sobre la secuencia original
        detalle_motivos_full[motivo] = pos_list

    # Construir cadena resumida para el TXT (mostrar hasta MAX_POS_SHOW por motivo)
    resumen_motivos = []
    for motivo, pos_list in detalle_motivos_full.items():
        if pos_list:
            if len(pos_list) > MAX_POS_SHOW:
                mostradas = pos_list[:MAX_POS_SHOW]
                resumen = f"{motivo}(@{','.join(map(str, mostradas))},...+{len(pos_list)-MAX_POS_SHOW})"
            else:
                resumen = f"{motivo}(@{','.join(map(str, pos_list))})"
        else:
            resumen = f"{motivo}()"
        resumen_motivos.append(resumen)

    resultados.append({
        "id": seq_id,
        "longitud": longitud,
        "traduccion": traduccion,
        # guardamos detalle completo como texto (posiciones separadas por ;)
        "motivos_full": ";".join([f"{m}:{','.join(map(str, detalle_motivos_full[m]))}" for m in MOTIVOS]),
        "motivos_preview": "; ".join(resumen_motivos)
    })

# Escribir archivo TXT (legible) en UTF-8
txt_path = "resultados_sesion4.txt"
with open(txt_path, "w", encoding="utf-8") as txt:
    header = f"{'ID':40} {'Long':6} {'Traducción (inicio...)':35} Motivos\n"
    txt.write(header)
    txt.write("-" * 120 + "\n")
    print(header.strip())
    print("-" * 120)
    for r in resultados:
        traduccion_preview = (r["traduccion"][:32] + "...") if len(r["traduccion"]) > 35 else r["traduccion"]
        line = f"{r['id'][:40]:40} {r['longitud']:6} {traduccion_preview:35} {r['motivos_preview']}\n"
        txt.write(line)
        print(line.strip())

# Escribir CSV completo en UTF-8 (ideal para Excel / importación)
csv_path = "resultados_sesion4.csv"
csv_fields = ["id", "longitud", "traduccion", "motivos_full"]
with open(csv_path, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_fields)
    writer.writeheader()
    for r in resultados:
        writer.writerow({
            "id": r["id"],
            "longitud": r["longitud"],
            "traduccion": r["traduccion"],
            "motivos_full": r["motivos_full"]
        })

print(f"\nArchivos generados: {txt_path}   {csv_path}")
print(f"Estrategia usada para codón parcial: {STRATEGY}")
