"""
ab1Scope: Interactive Analyzer for Sanger Electropherograms (ab1 files)

This script provides tools for:
- Interactive visualization and filtering of chromatogram peaks
- Simulated base calling with adjustable intensity/dominance thresholds
- Export of filtered sequences to FASTA and CSV
- Quality assessment using Phred-like scores
- Detection of ambiguous base calls and potential indels
- Multi-sample alignment and polymorphism detection
- Automated batch export of IUPAC-encoded sequences

Designed for use in Jupyter Notebooks or as part of larger pipelines.
Developed by Ricardo M. Borges and collaborators.
"""

import os
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from IPython.display import display
from ipywidgets import (
    Dropdown, IntSlider, FloatSlider, Checkbox, IntRangeSlider,
    VBox, HBox, Label, Output, interactive_output
)

# ==== Interactive AB1 Processing Function ====
def process_ab1(ab1_dir=".", return_data=False):
    ab1_files = [f for f in os.listdir(ab1_dir) if f.lower().endswith(".ab1")]
    if not ab1_files:
        raise FileNotFoundError("No .ab1 files found in the specified directory.")

    out = Output()

    # Store last run data in outer scope
    state = {}

    file_selector = Dropdown(options=ab1_files, description="File:")
    int_min_slider = IntSlider(value=10, min=0, max=1000, step=10, description="Min. Intensity")
    int_max_slider = IntSlider(value=20000, min=500, max=60000, step=100, description="Max. Intensity")
    dom_ratio_slider = FloatSlider(value=1.0, min=1.0, max=5.0, step=0.1, description="Dominance")
    norm_check = Checkbox(value=False, description="Normalize")
    smooth_check = Checkbox(value=True, description="Smooth")
    smooth_window_slider = IntSlider(value=5, min=1, max=25, step=2, description="Window")
    zoom_slider = IntRangeSlider(value=[0, 1000], min=0, max=2000, step=100, description="Zoom")

    controls = VBox([
        file_selector,
        VBox([
            HBox([int_min_slider, Label("Minimum total intensity to accept a peak.")]),
            HBox([int_max_slider, Label("Maximum total intensity allowed (remove saturated peaks).")]),
            HBox([dom_ratio_slider, Label("Minimum dominance ratio to consider a clean peak.")]),
            HBox([norm_check, Label("Normalize all signals (0 to 1).")]),
            HBox([smooth_check, Label("Apply smoothing to signals.")]),
            HBox([smooth_window_slider, Label("Size of smoothing window.")]),
            HBox([zoom_slider, Label("Zoom range (start and end of X-axis).")])
        ])
    ])

    def run_analysis(
        file_name, intensity_min_threshold, intensity_max_threshold,
        dominance_ratio_threshold, normalize, smooth, smooth_window, zoom_range
    ):
        ab1_path = os.path.join(ab1_dir, file_name)
        record = SeqIO.read(ab1_path, "abi")
        raw = record.annotations['abif_raw']

        A = np.array(raw['DATA9'], dtype=float)
        C = np.array(raw['DATA10'], dtype=float)
        G = np.array(raw['DATA11'], dtype=float)
        T = np.array(raw['DATA12'], dtype=float)
        base_positions = np.array(raw['PLOC1'], dtype=int)

        def smooth_signal(signal, window):
            return np.convolve(signal, np.ones(window)/window, mode='same')

        def normalize_signal(signal):
            return (signal - np.min(signal)) / (np.max(signal) - np.min(signal))

        if smooth:
            A = smooth_signal(A, smooth_window)
            C = smooth_signal(C, smooth_window)
            G = smooth_signal(G, smooth_window)
            T = smooth_signal(T, smooth_window)
        if normalize:
            A = normalize_signal(A)
            C = normalize_signal(C)
            G = normalize_signal(G)
            T = normalize_signal(T)

        filtered_seq = ""
        accepted_positions = []
        accepted_bases = []
        rejected_positions = []

        print("\nFirst 10 diagnostics:\n")
        for i, pos in enumerate(base_positions):
            signals = np.array([A[pos], C[pos], G[pos], T[pos]])
            total = np.sum(signals)
            top2 = np.sort(signals)[-2:]
            ratio = top2[0] / max(top2[1], 1e-6)
            base = "ACGT"[np.argmax(signals)]

            if i < 10:
                print(f"Position {pos} | signals={signals.round(2)} | total={total:.2f} | ratio={ratio:.2f} | base={base}")

            if total < intensity_min_threshold or total > intensity_max_threshold:
                filtered_seq += "N"
                rejected_positions.append(pos)
            else:
                filtered_seq += base
                accepted_positions.append(pos)
                accepted_bases.append(base)

        x = np.arange(len(A))
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x, y=A, mode='lines', name='A', line=dict(color='green')))
        fig.add_trace(go.Scatter(x=x, y=C, mode='lines', name='C', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=x, y=G, mode='lines', name='G', line=dict(color='black')))
        fig.add_trace(go.Scatter(x=x, y=T, mode='lines', name='T', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=accepted_positions, y=[1]*len(accepted_positions),
                                 mode='markers', name='Accepted', marker=dict(color='lime', size=6)))
        fig.add_trace(go.Scatter(x=rejected_positions, y=[1]*len(rejected_positions),
                                 mode='markers', name='Rejected', marker=dict(color='red', size=4)))

        if accepted_positions:
            y_text = [max(A[p], C[p], G[p], T[p]) + 10 for p in accepted_positions]
            fig.add_trace(go.Scatter(
                x=accepted_positions,
                y=y_text,
                text=accepted_bases,
                mode="text",
                textfont=dict(size=16, family="Courier New"),
                showlegend=False
            ))

        fig.update_layout(
            title=f"Electropherogram with Base Calls - {file_name}",
            xaxis_title="Data Point",
            yaxis_title="Intensity",
            xaxis_range=zoom_range
        )

        out.clear_output()
        with out:
            display(fig)
            print("\nFiltered Sequence (first 200 bases):")
            print(filtered_seq[:100])

            save_dir = os.path.join(ab1_dir, "data")
            os.makedirs(save_dir, exist_ok=True)
            csv_path = os.path.join(save_dir, file_name.replace(".ab1", "_sequence.csv"))

            df_sequence = pd.DataFrame({
                "position": accepted_positions,
                "base": accepted_bases,
                "A": [A[p] for p in accepted_positions],
                "C": [C[p] for p in accepted_positions],
                "G": [G[p] for p in accepted_positions],
                "T": [T[p] for p in accepted_positions],
                "total_intensity": [A[p] + C[p] + G[p] + T[p] for p in accepted_positions],
                "dominance_ratio": [
                    max(A[p], C[p], G[p], T[p]) / (sorted([A[p], C[p], G[p], T[p]])[-2] + 1e-6)
                    for p in accepted_positions
                ]
            })
            df_sequence.to_csv(csv_path, index=False)
            print(f"\nSequence CSV saved: {csv_path}")

            # Export FASTA
            fasta_path = os.path.join(save_dir, file_name.replace(".ab1", ".fasta"))
            with open(fasta_path, "w") as fasta_out:
                fasta_out.write(f">{file_name.replace('.ab1', '')}\n")
                for i in range(0, len(filtered_seq), 70):
                    fasta_out.write(filtered_seq[i:i+70] + "\n")
            print(f"✅ FASTA file saved: {fasta_path}")

            # Store result in state dict
            state['result'] = (base_positions, A, C, G, T)
            state['filtered_seq'] = filtered_seq
            state[file_name] = {
                "params": {
                    "intensity_min_threshold": intensity_min_threshold,
                    "intensity_max_threshold": intensity_max_threshold,
                    "dominance_ratio_threshold": dominance_ratio_threshold,
                    "normalize": normalize,
                    "smooth": smooth,
                    "smooth_window": smooth_window,
                    "zoom_range": zoom_range
                },
                "result": (base_positions, A, C, G, T),
                "filtered_seq": filtered_seq
            }

    widget_dict = {
        'file_name': file_selector,
        'intensity_min_threshold': int_min_slider,
        'intensity_max_threshold': int_max_slider,
        'dominance_ratio_threshold': dom_ratio_slider,
        'normalize': norm_check,
        'smooth': smooth_check,
        'smooth_window': smooth_window_slider,
        'zoom_range': zoom_slider
    }

    interactive_plot = interactive_output(run_analysis, widget_dict)
    display(controls, interactive_plot, out)

    if return_data:
        return state




# ==== Plot Quality Analysis ====
def plot_quality_analysis(base_positions, A, C, G, T, intensity_min_threshold=100, intensity_max_threshold=30000, file_name="output"):
    os.makedirs("images", exist_ok=True)
    dominance_ratios = []
    total_intensities = []
    phred_scores = []
    for pos in base_positions:
        signals = [A[pos], C[pos], G[pos], T[pos]]
        total = sum(signals)
        sorted_signals = sorted(signals, reverse=True)
        ratio = sorted_signals[0] / (sorted_signals[1] + 1e-6)
        phred = min(99, max(0, -10 * np.log10(1.0 / (ratio + 1e-6))))
        dominance_ratios.append(ratio)
        total_intensities.append(total)
        phred_scores.append(phred)
    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=base_positions, y=dominance_ratios, mode='lines', name='Dominance Ratio', line=dict(color='purple')))
    fig1.add_shape(type="line", x0=min(base_positions), x1=max(base_positions), y0=1.5, y1=1.5, line=dict(color="gray", dash="dash"))
    fig1.update_layout(title="Dominance Ratio (Top1 / Top2)", xaxis_title="Base position", yaxis_title="Ratio")
    fig1.write_html(f"images/{file_name}_dominance_ratio.html")
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=base_positions, y=total_intensities, mode='lines', name='Total Intensity', line=dict(color='black')))
    fig2.add_shape(type="rect", x0=min(base_positions), x1=max(base_positions), y0=0, y1=intensity_min_threshold, fillcolor='gray', opacity=0.2, line_width=0)
    fig2.add_shape(type="rect", x0=min(base_positions), x1=max(base_positions), y0=intensity_max_threshold, y1=max(total_intensities)+1000, fillcolor='orange', opacity=0.2, line_width=0)
    fig2.update_layout(title="Total Signal Intensity", xaxis_title="Base position", yaxis_title="Intensity")
    fig2.write_html(f"images/{file_name}_total_intensity.html")
    fig3 = go.Figure()
    fig3.add_trace(go.Scatter(x=base_positions, y=phred_scores, mode='lines', name='Quality Score', line=dict(color='teal')))
    fig3.update_layout(title="Phred-like Quality Score", xaxis_title="Base position", yaxis_title="Score (0–99)")
    fig3.write_html(f"images/{file_name}_quality_score.html")
    print(f"✅ HTMLs saved in ./images:\n- {file_name}_dominance_ratio.html\n- {file_name}_total_intensity.html\n- {file_name}_quality_score.html")

# ==== Analyze Quality Score Phred-like ====
def analyze_quality_score_phred_like(base_positions, A, C, G, T, file_name, ab1_dir, min_quality_threshold=15, min_stretch_length=5):
    signals = np.vstack([A[base_positions], C[base_positions], G[base_positions], T[base_positions]]).T
    dominance_ratios = np.array([max(s) / (sorted(s)[-2] + 1e-6) for s in signals])
    phred_scores = np.clip(-10 * np.log10(1 / (dominance_ratios + 1e-6)), 0, 60)
    low_quality_mask = phred_scores < min_quality_threshold
    low_regions = []
    start = None
    for i, low in enumerate(low_quality_mask):
        if low:
            if start is None:
                start = i
        elif start is not None:
            if i - start >= min_stretch_length:
                low_regions.append((start, i - 1))
            start = None
    if start is not None and len(base_positions) - start >= min_stretch_length:
        low_regions.append((start, len(base_positions) - 1))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=base_positions, y=phred_scores, mode='lines+markers', name='Phred-like Q', line=dict(color='black'), marker=dict(size=6)))
    for (start_idx, end_idx) in low_regions:
        fig.add_vrect(x0=base_positions[start_idx], x1=base_positions[end_idx], fillcolor='gray', opacity=0.3, line_width=0, annotation_text="Low Quality", annotation_position="top left")
    fig.update_layout(title=f"Simulated Phred-like Quality Scores - {file_name}", xaxis_title="Base Position", yaxis_title="Quality Score (Phred-like)", yaxis=dict(range=[0, 60]))
    output_dir = os.path.join(ab1_dir, "images")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, file_name.replace(".ab1", "_quality.html"))
    fig.write_html(output_path)
    print(f"✅ Quality analysis plot saved to: {output_path}")
    return phred_scores, dominance_ratios, low_regions

# ==== Detect Ambiguities and Indels ====
def detect_ambiguities_and_indels(base_positions, A, C, G, T, iupac_seq, file_name, ab1_dir, jump_threshold=20):
    output = []
    bases_iupac = list(iupac_seq)
    base_labels = "ACGT"
    for i, pos in enumerate(base_positions):
        base = bases_iupac[i]
        reason = ""
        if base not in base_labels:
            reason = "Ambiguous base (IUPAC)"
        if i > 0 and (base_positions[i] - base_positions[i-1]) > jump_threshold:
            reason += " + Possible indel (large jump)"
        if reason:
            signal = [A[pos], C[pos], G[pos], T[pos]]
            output.append({
                "position": pos,
                "IUPAC_base": base,
                "A": signal[0],
                "C": signal[1],
                "G": signal[2],
                "T": signal[3],
                "total_intensity": sum(signal),
                "dominance_ratio": max(signal) / (sorted(signal)[-2] + 1e-6),
                "reason": reason.strip(" +")
            })
    df = pd.DataFrame(output)
    output_dir = os.path.join(ab1_dir, "data")
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, file_name.replace(".ab1", "_ambiguities_and_indels.csv"))
    df.to_csv(csv_path, index=False)
    print(f"✅ Ambiguities and indels exported to: {csv_path}")
    return df

# ==== Export IUPAC FASTA ====
def export_iupac_fasta(iupac_seq, file_name, ab1_dir="."):
    os.makedirs(os.path.join(ab1_dir, "data"), exist_ok=True)
    record = SeqRecord(Seq(iupac_seq), id=file_name.replace(".ab1", ""), description="IUPAC Base Calls")
    fasta_path = os.path.join(ab1_dir, "data", file_name.replace(".ab1", "_iupac.fasta"))
    with open(fasta_path, "w") as f:
        SeqIO.write(record, f, "fasta")
    print(f"✅ IUPAC FASTA saved: {fasta_path}")

# ==== Export FASTA and Diagnostics CSV ====
def export_fasta_and_diagnostics_csv(file_name, ab1_dir, base_positions, A, C, G, T, accepted_positions, accepted_bases, filtered_seq, intensity_min_threshold, intensity_max_threshold):
    os.makedirs(os.path.join(ab1_dir, "data"), exist_ok=True)
    fasta_seq = "".join(accepted_bases)
    record = SeqRecord(Seq(fasta_seq), id=file_name.replace(".ab1", ""), description="Accepted base calls")
    fasta_path = os.path.join(ab1_dir, "data", file_name.replace(".ab1", "_accepted.fasta"))
    with open(fasta_path, "w") as f:
        SeqIO.write(record, f, "fasta")
    print(f"✅ FASTA saved: {fasta_path}")

    all_data = []
    for i, pos in enumerate(base_positions):
        signal_vec = [A[pos], C[pos], G[pos], T[pos]]
        total = sum(signal_vec)
        sorted_vals = sorted(signal_vec)
        ratio = sorted_vals[-1] / (sorted_vals[-2] + 1e-6)
        max_base = "ACGT"[np.argmax(signal_vec)]

        if total < intensity_min_threshold:
            reason = "low_total"
        elif total > intensity_max_threshold:
            reason = "saturated"
        elif ratio < 1.5:
            reason = "low_dominance"
        else:
            reason = "accepted"

        all_data.append({
            "position": pos,
            "A": signal_vec[0],
            "C": signal_vec[1],
            "G": signal_vec[2],
            "T": signal_vec[3],
            "total_intensity": total,
            "dominance_ratio": ratio,
            "base_call": max_base if reason == "accepted" else "N",
            "filter_reason": reason
        })

    df_diag = pd.DataFrame(all_data)
    csv_path = os.path.join(ab1_dir, "data", file_name.replace(".ab1", "_diagnostics.csv"))
    df_diag.to_csv(csv_path, index=False)
    print(f"✅ Diagnostics CSV saved: {csv_path}")

# ==== Compare AB1 Samples ====
def compare_ab1_samples(ab1_dir=".", file_list=None, intensity_threshold=10, max_threshold=20000, smooth=True, smooth_window=5):
    if file_list is None:
        file_list = [f for f in os.listdir(ab1_dir) if f.lower().endswith(".ab1")]
    if len(file_list) < 2:
        print("!! Need at least two .ab1 files for comparison.")
        return

    sequences_dict = {}
    position_sets = {}
    trace_dict = {}

    def smooth_signal(signal, window):
        return np.convolve(signal, np.ones(window)/window, mode='same')

    def extract_filtered_sequence(file_name):
        record = SeqIO.read(os.path.join(ab1_dir, file_name), "abi")
        raw = record.annotations['abif_raw']
        A = np.array(raw['DATA9'], dtype=float)
        C = np.array(raw['DATA10'], dtype=float)
        G = np.array(raw['DATA11'], dtype=float)
        T = np.array(raw['DATA12'], dtype=float)
        base_positions = np.array(raw['PLOC1'], dtype=int)
        if smooth:
            A, C, G, T = [smooth_signal(s, smooth_window) for s in (A, C, G, T)]
        seq_dict = {}
        for i, pos in enumerate(base_positions):
            signals = np.array([A[pos], C[pos], G[pos], T[pos]])
            total = np.sum(signals)
            if total < intensity_threshold or total > max_threshold:
                seq_dict[pos] = "N"
            else:
                base = "ACGT"[np.argmax(signals)]
                seq_dict[pos] = base
        return seq_dict, A, C, G, T

    fig = go.Figure()
    common_positions = None
    for fname in file_list:
        seq_dict, A, C, G, T = extract_filtered_sequence(fname)
        sequences_dict[fname] = seq_dict
        position_sets[fname] = set(seq_dict.keys())
        trace_dict[fname] = {"A": A, "C": C, "G": G, "T": T}

        x = np.arange(len(A))
        for nt, signal in zip("ACGT", (A, C, G, T)):
            fig.add_trace(go.Scatter(x=x, y=signal, mode='lines', name=f"{fname} - {nt}", line=dict(width=1)))

        common_positions = position_sets[fname] if common_positions is None else common_positions & position_sets[fname]

    if not common_positions:
        print("!! No common base positions to compare across samples.")
        return

    sorted_common_positions = sorted(common_positions)
    aligned_df = pd.DataFrame({pos: {fname: sequences_dict[fname][pos] for fname in file_list} for pos in sorted_common_positions}).T
    aligned_df.index.name = "Base_Position"
    aligned_df["Polymorphism"] = aligned_df.nunique(axis=1).gt(1)

    print("\nFirst 20 polymorphic positions:")
    display(aligned_df[aligned_df["Polymorphism"]].head(20))

    os.makedirs(os.path.join(ab1_dir, "images"), exist_ok=True)
    fig.write_html(os.path.join(ab1_dir, "images", "comparison_plot.html"))
    print("✅ Comparison plot saved as HTML.")

    aligned_df.to_csv(os.path.join(ab1_dir, "data", "sequence_alignment.csv"))
    print("✅ Alignment saved to 'data/sequence_alignment.csv'")

# ==== Generate and Export IUPAC for Files ====
def call_iupac_sequence(base_positions, A, C, G, T, intensity_threshold=10, max_threshold=60000):
    iupac_codes = {
        (True, False, False, False): 'A', (False, True, False, False): 'C',
        (False, False, True, False): 'G', (False, False, False, True): 'T',
        (True, True, False, False): 'M', (True, False, True, False): 'R',
        (True, False, False, True): 'W', (False, True, True, False): 'S',
        (False, True, False, True): 'Y', (False, False, True, True): 'K',
        (True, True, True, False): 'V', (True, True, False, True): 'H',
        (True, False, True, True): 'D', (False, True, True, True): 'B',
        (True, True, True, True): 'N',
    }
    sequence = ""
    for pos in base_positions:
        signals = np.array([A[pos], C[pos], G[pos], T[pos]])
        total = np.sum(signals)
        if total < intensity_threshold or total > max_threshold:
            sequence += 'N'
        else:
            top_value = max(signals)
            presence = tuple((signals / (top_value + 1e-6)) > 0.33)
            sequence += iupac_codes.get(presence, 'N')
    return sequence

def generate_and_export_iupac_for_files(file_list, ab1_dir, state=None):
    """
    Gera sequência com código IUPAC e exporta para FASTA para todos os arquivos .ab1.
    Pode opcionalmente usar o `state` com dados já carregados.
    """
    if state is None:
        state = {}

    for fname in file_list:
        ab1_path = os.path.join(ab1_dir, fname)
        print(f"Processando: {fname}")

        # Se já estiver no state, reaproveita
        if fname in state and "result" in state[fname]:
            base_pos, A, C, G, T = state[fname]["result"]
        else:
            record = SeqIO.read(ab1_path, "abi")
            raw = record.annotations['abif_raw']
            A = np.array(raw['DATA9'], dtype=float)
            C = np.array(raw['DATA10'], dtype=float)
            G = np.array(raw['DATA11'], dtype=float)
            T = np.array(raw['DATA12'], dtype=float)
            base_pos = np.array(raw['PLOC1'], dtype=int)

        # Gera sequência IUPAC
        iupac_seq = call_iupac_sequence(base_pos, A, C, G, T)

        # Exporta para FASTA
        export_iupac_fasta(iupac_seq, fname, ab1_dir)

        # Atualiza o state (opcional)
        state[fname] = state.get(fname, {})
        state[fname]["result"] = (base_pos, A, C, G, T)
        state[fname]["iupac_seq"] = iupac_seq

    print("✅ Todas as sequências foram exportadas para FASTA com código IUPAC.")
    return state
