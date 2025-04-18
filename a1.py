import streamlit as st
import numpy as np
import os

# Constants
g = 9.8
pi = 3.14
rho = 1.025

def run_calculations(cc, ratios, L, draft, data_folder='data', damping_folder='DAMPING'):
    n = len(cc)
    m = len(cc[0])
    wsq = g * 2 * pi / L / 1.2

    files = os.listdir(data_folder)
    ratio_in = {8: 7, 10: 7}
    for idx, i in enumerate(ratios):
        ratio_in[i * 2.5] = idx

    sa, beta, b_by_t, xaxis = [], [], [], []
    a_cf, d_cf, wpa = [], [], []

    for i in range(1, m):
        pr = pr2 = ar = b = 0
        for j in range(1, n):
            ar += (cc[j][i] + pr) * (cc[j][0] - pr2)
            pr, pr2 = cc[j][i], cc[j][0]
            b = max(b, cc[j][i])
        sa.append(ar)
        beta_val = round(ar / (b * 2 * draft), 1)
        beta.append(max(beta_val, 0.5))
        b_by_t_val = b * 2 / draft
        b_by_t.append(b_by_t_val)
        xaxis_val = wsq * b / g
        xaxis.append(xaxis_val)

        beta_key = int(beta[-1] * 10)
        ratio_key = max(1, round(b_by_t_val * 2.5, 0))
        tem = ratio_in.get(ratio_key, 0)

        # Added mass coefficient
        with open(f'{data_folder}/{beta_key}.txt') as f:
            con = [[0 if j == '' else float(j) for j in i.split('\t')] for i in f.read().split('\n')]
        dd, y = 1000, -1
        for row in con:
            if abs(xaxis_val - row[tem * 2]) < dd:
                dd, y = abs(xaxis_val - row[tem * 2]), row[tem * 2 + 1]
        a_cf.append(y)

        # Damping coefficient
        with open(f'{damping_folder}/{beta_key}.txt') as f:
            con = [[0 if j == '' else float(j) for j in i.split('\t')] for i in f.read().split('\n')]
        dd, y = 1000, -1
        for row in con:
            if abs(xaxis_val - row[tem * 2]) < dd:
                dd, y = abs(xaxis_val - row[tem * 2]), row[tem * 2 + 1]
        d_cf.append(y)

    aa = list(a_cf[::-1])
    aa.append(0)
    aa = list(aa[::-1])
    a3 = a55 = i_mass = 0
    for idx, i in enumerate(cc[0][1:]):
        delta = i - cc[0][idx]
        a3 += (aa[idx] + aa[idx + 1]) / 2 * delta * sa[idx] * rho
        a55 += (aa[idx] + aa[idx + 1]) / 2 * delta * sa[idx] * rho * (L / 2 - i) ** 2
        i_mass += delta * sa[idx] * rho * (L / 2 - i) ** 2

    aa = list(d_cf[::-1])
    aa.append(0)
    aa = list(aa[::-1])
    b3 = b55 = 0
    for idx, i in enumerate(cc[0][1:]):
        delta = i - cc[0][idx]
        b3 += (aa[idx] + aa[idx + 1]) / 2 * delta * sa[idx] * rho
        b55 += (aa[idx] + aa[idx + 1]) / 2 * delta * sa[idx] * rho * (L / 2 - i) ** 2

    for i in range(1, n):
        pr = pr2 = ar = 0
        for j in range(1, m):
            ar += (cc[i][j] + pr) * (cc[0][j] - pr2)
            pr, pr2 = cc[i][j], cc[0][j]
        wpa.append(ar)

    return {
        "A33": a3,
        "A55": a55,
        "B33": b3,
        "B55": b55,
        "C33": wpa[0] * rho * g if wpa else 0,
        "C55": wpa[-1] * rho * g if wpa else 0
    }

# Streamlit App
st.set_page_config(page_title="Ship Hydrodynamics Calculator", layout="centered")
st.title("🚢 Ship Hydrodynamics Calculator")

st.sidebar.header("Ship Parameters")
L = st.sidebar.number_input("Enter Ship Length (L) in meters", value=122.0, step=0.1)
draft = st.sidebar.number_input("Enter Ship Draft in meters", value=7.8, step=0.1)

uploaded_file = st.file_uploader("Upload Offset Table (.txt)", type=["txt"])

if uploaded_file:
    content = uploaded_file.read().decode('utf-8')
    cc = [[float(j) for j in i.split('\t')] for i in content.strip().split('\n')]

    # Load ratios
    with open("data.txt", 'r') as f:
        ratios = [float(i) for i in f.read().split('\t\t')]

    results = run_calculations(cc, ratios, L, draft)

    st.header("Results")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("A33 (Added mass for heave)", f"{results['A33']:.3f}")
        st.metric("B33 (Damping for heave)", f"{results['B33']:.3f}")
        st.metric("C33 (Restoring force for heave)", f"{results['C33']:.3f}")
    with col2:
        st.metric("A55 (Added mass for pitch)", f"{results['A55']:.3f}")
        st.metric("B55 (Damping for pitch)", f"{results['B55']:.3f}")
        st.metric("C55 (Restoring force for pitch)", f"{results['C55']:.3f}")
else:
    st.info("Please upload an offset table in `.txt` format.")
