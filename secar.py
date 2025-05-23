import os

# Read input data from file
with open(f'in.txt', 'r') as f:
    cc = [[float(j) for j in i.split('\t')] for i in f.read().strip().split('\n')]

n = len(cc)  # Number of rows
m = len(cc[0])  # Number of columns

# Constants
g = 9.8  # Gravitational acceleration (m/s²)
pi = 3.14  # Pi approximation
draft = 7.8  # Ship draft (m)
L = 122  # Ship length (m)
rho = 1.025  # Water density (kg/m³), used in later calculations

# Calculate wave frequency squared
wsq = g * 2 * pi / L / 1.2

# Read ratio data
with open('data.txt', 'r') as f:
    ratios = [float(i) for i in f.read().split('\t\t')]

# List available data files
files = os.listdir('data')

# Initialize ratio lookup dictionary
ratio_in = {8: 7, 10: 7}
for idx, i in enumerate(ratios):
    ratio_in[i * 2.5] = idx 

# Initialize arrays for calculations
sa = []      # Sectional areas
beta = []    # Shape coefficient
b_by_t = []  # Beam to draft ratio
xaxis = []   # Frequency parameter
wpa = []     # Waterplane areas
a_cf = []    # Added mass coefficients
d_cf = []    # Damping coefficients

# Process each column (ship section)
for i in range(1, m):
    pr = 0    # Previous value
    pr2 = 0   # Previous x-coordinate
    ar = 0    # Area
    b = 0     # Maximum breadth
    
    # Calculate sectional area and find maximum breadth
    for j in range(1, n):
        ar += (cc[j][i] + pr) * (cc[j][0] - pr2)  # Area calculation using trapezoidal rule
        pr = cc[j][i]   # Update previous y-value
        pr2 = cc[j][0]  # Update previous x-value
        b = max(b, cc[j][i])  # Update maximum breadth
    
    # Store calculated values
    sa.append(ar)  # Sectional area
    beta.append(round(sa[-1] / (b * 2 * draft), 1))  # Shape coefficient
    
    # Ensure minimum beta value
    if beta[-1] < 0.5:
        beta[-1] = 0.5
    
    b_by_t.append(b * 2 / draft)  # Beam to draft ratio
    xaxis.append(wsq * b / g)  # Frequency parameter
    
    # Read added mass coefficients from data files
    with open(f'data/{int(beta[-1] * 10)}.txt', 'r') as f:
        con = [[0 if j == '' else float(j) for j in i.split('\t')] for i in f.read().split('\n')]
        
        # Find appropriate ratio index
        tem = ratio_in[max(1, round(b_by_t[-1] * 2.5, 0))]
        
        # Find closest matching value
        dd = 1000  # Initialize with a large value
        y = -1     # Default value
        
        for row in con:
            if abs(xaxis[-1] - row[tem * 2]) < dd:
                dd = abs(xaxis[-1] - row[tem * 2])
                y = row[tem * 2 + 1]
        
        a_cf.append(y)  # Store added mass coefficient
    
    # Read damping coefficients from data files
    with open(f'DAMPING/{int(beta[-1] * 10)}.txt', 'r') as f:
        con = [[0 if j == '' else float(j) for j in i.split('\t')] for i in f.read().split('\n')]
        
        # Find appropriate ratio index
        tem = ratio_in[max(1, round(b_by_t[-1] * 2.5, 0))]
        
        # Find closest matching value
        dd = 1000  # Initialize with a large value
        y = -1     # Default value
        
        for row in con:
            if abs(xaxis[-1] - row[tem * 2]) < dd:
                dd = abs(xaxis[-1] - row[tem * 2])
                y = row[tem * 2 + 1]
        
        d_cf.append(y)  # Store damping coefficient


aa = list(a_cf[::-1])
aa.append(0)
aa = list(aa[::-1])

# Initialize variables for heave calculations
p1 = 0       # Unused variable
a3 = 0       # Heave added mass
a55 = 0      # Pitch added mass moment
i_mass = 0   # Mass moment of inertia

# Calculate heave and pitch parameters
for idx, i in enumerate(cc[0][1:]):
    a3 += (aa[idx] + aa[idx + 1]) / 2 * (i - cc[0][idx]) * sa[idx] * rho
    a55 += (aa[idx] + aa[idx + 1]) / 2 * (i - cc[0][idx]) * sa[idx] * rho * (L / 2 - i) ** 2
    i_mass += (i - cc[0][idx]) * sa[idx] * rho * (L / 2 - i) ** 2

aa = list(d_cf[::-1])
aa.append(0)
aa = list(aa[::-1])

# Initialize variables for damping calculations
p1 = 0  # Unused variable
b3 = 0  # Heave damping
b55 = 0 # Pitch damping moment

# Calculate damping coefficients
for idx, i in enumerate(cc[0][1:]):
    b3 += (aa[idx] + aa[idx + 1]) / 2 * (i - cc[0][idx]) * sa[idx] * rho
    b55 += (aa[idx] + aa[idx + 1]) / 2 * (i - cc[0][idx]) * sa[idx] * rho * (L / 2 - i) ** 2

# Calculate waterplane areas
for i in range(1, n):
    pr = 0    # Previous value
    ar = 0    # Area
    b = 0     # Maximum breadth (unused)
    pr2 = 0   # Previous x-coordinate
    
    for j in range(1, m):
        ar += (cc[i][j] + pr) * (cc[0][j] - pr2)  # Area calculation using trapezoidal rule
        pr = cc[i][j]   # Update previous y-value
        pr2 = cc[0][j]  # Update previous x-value
    
    wpa.append(ar)  # Store waterplane area

# Write results to output file
with open('output.txt', 'w') as f:
    print('SEC AREA\t' + '\t'.join([str(round(i, 2)) for i in sa]), file=f)
    print('Beta\t' + '\t'.join([str(round(i, 2)) for i in beta]), file=f)
    print('B/T\t' + '\t'.join([str(round(i, 2)) for i in b_by_t]), file=f)
    print('X_Axis\t' + '\t'.join([str(round(i, 2)) for i in xaxis]), file=f)
    print('Added mass coeff\t' + '\t'.join([str(round(i, 2)) for i in a_cf]), file=f)
    print('Damping Coeff\t' + '\t'.join([str(round(i, 2)) for i in d_cf]), file=f)

    print(f'\nA33\t{a3:.3g}\nA55\t{a55:.3g}\nB33\t{b3:.4g}\nB55\t{b55:.3g}\ni_mass\t{i_mass:.3g}\n', file=f)

    print('\nWPA\tC33\n' + '\n'.join([f"{str(round(i, 2))}\t{str(round(i*rho*g, 2))}" for i in wpa]), file=f)
