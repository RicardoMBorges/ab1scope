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


### Comparação com Métodos Quantitativos de Proporção de Cópias

Embora o ab1scope ofereça uma plataforma interativa robusta para análise de qualidade e extração de sequências a partir de arquivos .ab1, seu foco está em curadoria manual assistida, visualização e detecção qualitativa de sinais dominantes. Ele não realiza estimativas quantitativas precisas de proporções de variantes (como SNPs ou edições CRISPR), como proposto em metodologias avançadas descritas por Seroussi (2021).
O que o ab1scope oferece:

- Interface Jupyter interativa com ajustes em tempo real

- Diagnóstico visual e exportação de sequências filtradas (FASTA/CSV)

- Análise em lote com rastreamento de parâmetros

- Identificação de ambiguidade por dominância e codificação IUPAC

- Ideal para ensino, controle de qualidade, curadoria e pré-processamento

### O que não é escopo do ab1scope:

- Inferência quantitativa precisa de proporções de cópias

- Correção de distorções sistemáticas de fluorescência

- Integração com modelos de edição genômica (e.g., TIDE, ICE)

- Análise bidirecional combinada de arquivos forward/reverse

    Recomendação: Para aplicações como análise de eficiência de edição por CRISPR, detecção de heteroplasmia mitocondrial ou epigenética quantitativa, sugerimos ferramentas específicas como EditR, TIDE, ICE, BEAT, entre outras.

### 🛑 Limitações e Considerações

O ab1scope foi projetado para análise interativa e curadoria qualitativa de dados de sequenciamento Sanger. Embora seja uma ferramenta poderosa para inspeção visual, ensino e diagnóstico básico de qualidade, ele não substitui métodos quantitativos rigorosos em determinadas aplicações.


#### ❌ Limitações atuais:

- Não estima proporções alélicas ou variantes com precisão quantitativa.
Ideal para detectar ambiguidade e dominância, mas não calcula frações de alelos (ex.: 70% A, 30% G) como modelos estatísticos fazem.

- Não corrige efeitos de fluorescência ou sobreposição espectral.
Os sinais são tratados diretamente a partir das intensidades brutas ou suavizadas, sem deconvolução de canais.

- Não integra arquivos forward/reverse nem realiza montagem.
    Cada .ab1 é processado individualmente; não há sobreposição de leituras ou geração de consenso multiarquivo.

- Não detecta ou classifica edições CRISPR, SNPs, ou inserções/deleções com precisão de bioinformática.
    Há uma função básica de detecção de ambiguidades e indels visuais, mas sem inferência estatística de mutações.

- Não substitui ferramentas específicas como TIDE, ICE, EditR, BEAT.
