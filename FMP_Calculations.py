import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import math

#Constants
AVOGADRO_NUMBER = 6.022e23

def create_materials_table():
    conn = sqlite3.connect('material_database.db')
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS materials (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            material_name TEXT NOT NULL,
            element_composition TEXT NOT NULL,
            atomic_numbers TEXT NOT NULL,
            valence_electrons TEXT NOT NULL,
            molar_mass REAL NOT NULL,
            density REAL NOT NULL,
            band_gap REAL NOT NULL
        )
    ''')
    conn.commit()
    conn.close()

# Function to add materials into Database
def add_material(material_name, element_composition, atomic_numbers, valence_electrons, molar_mass, density, band_gap):
    conn = sqlite3.connect('material_database.db')
    cursor = conn.cursor()
    create_materials_table()
    # Check if the material already exists in the database
    cursor.execute('SELECT * FROM materials WHERE material_name = ?', (material_name,))
    existing_material = cursor.fetchone()
    
    if existing_material:
        print(f'Material "{material_name}" already exists in the database.')
    else:
        cursor.execute('''
            INSERT INTO materials (material_name, element_composition, atomic_numbers, valence_electrons, molar_mass, density, band_gap)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (material_name, element_composition, atomic_numbers, valence_electrons, molar_mass, density, band_gap))
        
        # Save changes in the database
        conn.commit()
        conn.close()
        print(f'Material "{material_name}" has been added to the database.')

# Example to Add Material
#add_material("MoS2", "Mo=1,S=2", "Mo=42,S=18","Mo=6,S=6", 160.07, 5.06, 1.23)

# Function to get material data from Database
def get_material_data(material_name):
    conn = sqlite3.connect('material_database.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM materials WHERE material_name = ?', (material_name,))
   # material_name, element_composition, atomic_numbers, molar_mass, density, band_gap 
    material_info = cursor.fetchone()
    conn.close()
    if material_info:
        return material_info
    else:
        print(f'No material called "{material_name}" in database')

def list_all_tables():
    conn = sqlite3.connect('material_database.db')
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()   
    conn.close()
    return [table[0] for table in tables]

def get_column_names(table):
    conn = sqlite3.connect('material_database.db')
    cursor = conn.cursor()
    
    cursor.execute(f"PRAGMA table_info({table})")
    columns = cursor.fetchall()
    conn.close() 
    if columns:
        column_names = [column[1] for column in columns]
        return column_names
    else:
        print(f'get_column_names - No table called "{table}" found in the database.')

#print(get_column_names("materials"))   

def get_all_material_names(table):
    conn = sqlite3.connect('material_database.db')
    cursor = conn.cursor()
    
    cursor.execute(f"SELECT material_name FROM {table}")
    materials = cursor.fetchall()  
    conn.close()
    if materials:
        material_names = [material[0] for material in materials]
        return material_names
    else:
        print(f'get_all_material_names - No table called "{table}" found in the database')


def fmp_seah_S1_with_material(material_name, E):
    id,material_name, element_composition, atomic_numbers, valence_electrons, molar_mass, density, band_gap = get_material_data(material_name)
    if element_composition is not None:
        sum_Z = 0
        sum_atoms = 0
        # Number of Atoms in Material
        for elements in element_composition.split(','):
            atom, N_atom = elements.split("=")
            sum_atoms += int(N_atom)
        # Claculate mean atomic number
        for element_info in atomic_numbers.split(','):
            element, Z = element_info.split('=')
            Z = int(Z)  # Convert atomic number (Z) to an integer
            sum_Z += Z

        Z = sum_Z / sum_atoms
        # calculate the free mean path lambda (l)
        a = (molar_mass * 1e21 / (density * sum_atoms * AVOGADRO_NUMBER))**(1/3)
        l = (4 + 0.44 * Z**0.5 + 0.104 * E**0.872) * a**1.7 / (Z**0.3 * (1 - 0.02 * band_gap))
        return l
    else:
        print(f'Material {material_name} not in Database.')
        return None

def tpp_2mfmp_with_material(material_name, E):
    id,material_name, element_composition, atomic_numbers, valence_electrons, molar_mass, density, band_gap = get_material_data(material_name)
    # Calculation of valence electrons
    N_v = 0
    for electron in valence_electrons.split(','):
        atom, N_electron = electron.split("=")
        N_v += int(N_electron)

    U = N_v * density / molar_mass
    E_p = 28.8 * math.sqrt(U)
    gamma = 0.191 * density**-0.50
    beta = -0.1 + 0.944 / math.sqrt(E_p**2 + band_gap**2) + 0.069 * density**0.1
    C = 1.97 - 0.91 * U
    D = 53.4 - 20.8 * U
    l = E / (E_p**2 * (beta * math.log(gamma * E+ 1e-10) - C / E + D / E**2)) # Added a small constant to avoid log(0)
    l_nm = l * 0.1
    return l_nm

def plot_results(material_name, start, end, points):
    fmp_results = []
    tpp_results = []
    E_values = np.linspace(start, end, points)  # You can adjust the number of points
    for E in E_values:
        fmp_result = fmp_seah_S1_with_material(material_name, E)
        tpp_result = tpp_2mfmp_with_material(material_name, E)

        # Exclude None results from the plot
        if fmp_result is not None:
            fmp_results.append(fmp_result)
        else:
            fmp_results.append(np.nan)  # Use np.nan to leave a gap in the plot
        
        if tpp_result is not None:
            tpp_results.append(tpp_result)
        else:
            tpp_results.append(np.nan)

    plt.figure(figsize=(10, 6))

    plt.plot(E_values, fmp_results, label='S1 Results')
    plt.plot(E_values, tpp_results, label='2M-TPP Results')

    plt.yscale('log')  # Set the y-axis to be logarithmic
    plt.xscale('log')  # Set the y-axis to be logarithmic
    plt.xlabel('Kinetic Energy [eV]')
    plt.ylabel('Inelastic mean free path [nm]')
    plt.title(f'Results for {material_name}')
    plt.legend()
    plt.grid(True)
    plt.show()

#plot_results("MoS2",60,2000,200)

# def main():
#     #No idea at the moment how to set up main()
#     print("Welcome \n")
#     conn = sqlite3.connect('material_database.db')
#     create_materials_table()
#     column_names = get_column_names()
#     material_names = get_all_material_names()
#     print("The materials table in the DataBase contains the following columns: " + ', '.join(column_names) + "\n")
#     print("Currently, the materials table contains the following materials and properties needed for calculations: " + ', '.join(material_names) + "\n")
# if __name__ == '__main__':
#     main()