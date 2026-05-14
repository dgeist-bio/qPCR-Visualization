import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path = 'qPCR_data.xlsx' # Leg die Datei einfach in den gleichen Ordner

try:
    df = pd.read_excel(path, skiprows=2)
    
    # 2. Daten extrahieren (Namen, ddCT, SD)
    df_plot = df.iloc[:, [0, 13, 4]].dropna() 
    df_plot.columns = ['Sample Name', 'ddCT', 'SD']
    
    # Numerische Konvertierung
    df_plot['ddCT'] = pd.to_numeric(df_plot['ddCT'], errors='coerce')
    df_plot['SD'] = pd.to_numeric(df_plot['SD'], errors='coerce')
    df_plot = df_plot.dropna()

    # 3. Viridis-Farben generieren
    # Erstellt einen Farbverlauf passend zur Anzahl der Balken
    colors = plt.cm.viridis(np.linspace(0, 0.8, len(df_plot)))

    # 4. Plotten
    plt.figure(figsize=(12, 7))
    
    bars = plt.bar(df_plot['Sample Name'], df_plot['ddCT'], 
                   yerr=df_plot['SD'], 
                   capsize=7, 
                   color=colors, # Hier ist dein Viridis!
                   edgecolor='black',
                   error_kw={'ecolor': 'black', 'alpha': 0.7, 'linewidth': 2})

    plt.axhline(1.0, color='red', linestyle='--', label='Baseline (DMSO)')

    plt.title('Relative Genexpression (qPCR Analysis)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Relative Expression ($\Delta\Delta C_q$)')
    plt.legend()
    plt.tight_layout()

    plt.show()
    
except FileNotFoundError:
    print("Bitte lege die Excel-Datei 'qPCR_data.xlsx' in den Ordner.")