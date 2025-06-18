# ab1scope - Visualização, Análise e Exportação de Dados de Sequenciamento Sanger (.ab1)

**An interactive toolkit for visualizing, analyzing, and exporting Sanger sequencing (.ab1) files in Jupyter Notebooks.**
(P) ab1scope é uma ferramenta interativa feita em Python (versão 3.11) para analisar arquivos de sequenciamento Sanger no formato .ab1. Ela permite visualizar eletroferogramas, aplicar critérios de qualidade, gerar sequências com códigos IUPAC, detectar ambiguidades e exportar resultados em .csv e .fasta.

> Built with Python 3.11 and designed for educational and analytical purposes in molecular biology, biotechnology, and bioinformatics labs.

---

## 📦 Features

- Interactive electropherogram viewer
- Adjustable filters: signal intensity, dominance ratio, smoothing, zoom
- Quality score estimation (Phred-like)
- Detection of ambiguous bases and potential indels
- IUPAC sequence generation and FASTA export
- Side-by-side comparison of multiple `.ab1` files
- Export of base calls and metadata to CSV

---
Como Usar
- Baixe o projeto completo (incluindo run.bat).
- Coloque seus arquivos .ab1 na pasta desejada.
- Dê dois cliques no arquivo run.bat para abrir automaticamente o Jupyter Notebook no navegador.

 No Notebook:

- Defina o caminho da pasta onde estão os arquivos .ab1
- Execute as células sequencialmente para visualizar, ajustar parâmetros e exportar os dados.
---

## 📁 Directory Structure
main/
│
├── your_ab1_files_directory/
│ ├── sample1.ab1
│ ├── sample2.ab1
│ ├── ...
| └── ab1scope_notebook.ipynb # Interactive Jupyter interface
├── ab1scope.py # Core module
└── README.md

