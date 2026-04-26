import argparse
import pandas as pd
import sys

def merge_condition(row, ptm_col, cond_col):
    """Helper to cleanly merge a PTM and Condition column if they exist."""
    ptm = str(row.get(ptm_col, '')).strip()
    cond = str(row.get(cond_col, '')).strip()
    
    if ptm and cond:
        return f"{ptm}-{cond}"
    elif ptm:
        return ptm
    else:
        return cond

def main():
    # 1. Set up Command Line Arguments
    parser = argparse.ArgumentParser(description="Generate a Bingo Matrix and single-variable comparisons.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input experiment.csv file")
    args = parser.parse_args()

    # 2. Load the Data
    try:
        df = pd.read_csv(args.input)
    except FileNotFoundError:
        print(f"Error: Could not find '{args.input}'. Please check the path.")
        sys.exit(1)

    # Clean up NaNs safely across all possible condition columns
    for col in ['ptm_a', 'condition_a', 'ptm_b', 'condition_b']:
        if col in df.columns:
            df[col] = df[col].fillna('')

    # 3. Create Merged Conditions for BOTH Protein A and Protein B
    df['merged_cond_a'] = df.apply(lambda x: merge_condition(x, 'ptm_a', 'condition_a'), axis=1)
    df['merged_cond_b'] = df.apply(lambda x: merge_condition(x, 'ptm_b', 'condition_b'), axis=1)

    # Define the Groups for the Bingo Sheet
    df['Group_A'] = df['chain_a'] + '_' + df['merged_cond_a']
    df['Group_B'] = df['chain_b'] + '_' + df['merged_cond_b']

    # Full label for the final output
    df['Full_Experiment'] = df['Group_A'] + ' + ' + df['Group_B']

    # ==========================================
    # OUTPUT 1: Create the "Bingo Matrix"
    # ==========================================
    bingo_matrix = pd.crosstab(df['Group_A'], df['Group_B']).replace({1: '✅', 0: '-'})

    print("\n" + "="*80)
    print("--- BINGO MATRIX ---")
    print("="*80)
    print(bingo_matrix.to_markdown())
    print("\n")

    bingo_matrix.to_csv('bingo_matrix.csv')

    # ==========================================
    # OUTPUT 2: Generate Statistical Comparisons
    # ==========================================
    # The biological state is now defined by perfectly symmetrical 4 pillars:
    base_conditions = ['chain_a', 'merged_cond_a', 'chain_b', 'merged_cond_b']
    
    # Friendly names for the output table
    display_names = {
        'chain_a': 'Protein A',
        'merged_cond_a': 'Condition A',
        'chain_b': 'Protein B',
        'merged_cond_b': 'Condition B'
    }
    
    comparisons = []
    experiments = df.to_dict('records')

    # Iterate through every unique pair of experiments
    for i in range(len(experiments)):
        for j in range(i + 1, len(experiments)):
            exp1 = experiments[i]
            exp2 = experiments[j]
            
            # Check how many of our 4 pillars differ
            differences = []
            for col in base_conditions:
                if exp1[col] != exp2[col]:
                    differences.append(col)
            
            # If EXACTLY ONE pillar changes, it's a valid 1:1 comparison
            if len(differences) == 1:
                changed_var = differences[0]
                comparisons.append({
                    'Experiment 1': exp1['Full_Experiment'],
                    'Experiment 2': exp2['Full_Experiment'],
                    'Changed Variable': display_names[changed_var],
                    'Exp 1 State': exp1[changed_var],
                    'Exp 2 State': exp2[changed_var]
                })

    # Convert list of valid comparisons into a DataFrame
    comparisons_df = pd.DataFrame(comparisons)

    print("="*80)
    print("--- VALID STATISTICAL COMPARISONS (1 Variable Changed) ---")
    print("="*80)
    if not comparisons_df.empty:
        print(comparisons_df.to_markdown(index=False))
        comparisons_df.to_csv('statistical_comparisons.csv', index=False)
        print("\n✅ Success! Saved to 'bingo_matrix.csv' and 'statistical_comparisons.csv'\n")
    else:
        print("No valid single-variable comparisons found.\n")

if __name__ == "__main__":
    main()
