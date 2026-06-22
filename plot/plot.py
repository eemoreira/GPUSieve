import argparse
import pandas as pd
import matplotlib.pyplot as plt

# argumentos da CLI
parser = argparse.ArgumentParser(
    description="Plota gráfico de benchmark do sieve"
)

parser.add_argument(
    "arquivo",
    help="Arquivo CSV de entrada"
)

parser.add_argument(
    "-o",
    "--output",
    default="grafico.png",
    help="Nome do arquivo PNG de saída"
)

args = parser.parse_args()

df = pd.read_csv(args.arquivo)
plt.figure(figsize=(10, 6))

plt.plot(
    df["N"],
    df["Byte Sieve Time (ms)"],
    marker="o",
    label="Byte Sieve"
)

plt.plot(
    df["N"],
    df["Bit Sieve Time (ms)"],
    marker="o",
    label="Bit Sieve"
)

plt.xscale("log", base=2)

plt.xlabel("N")
plt.ylabel("Tempo (ms)")
plt.title("Byte Sieve vs Bit Sieve")

plt.grid(True, which="both", linestyle="--")
plt.legend()

plt.tight_layout()

plt.savefig(args.output, dpi=300)
plt.show()
