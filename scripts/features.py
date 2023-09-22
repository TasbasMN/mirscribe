def generate_alignment_string_from_dot_bracket(df):
    full_strings = []
    for _, row in df.iterrows():
        start_string = (row.mirna_start) * "0"
        mid_string = row["mirna_dot_bracket_5to3"].replace(".", "0").replace(")", "1")
        end_string = (len(row.mirna_sequence) - row.mirna_end -1) * "0"
        
        full_string = start_string + mid_string + end_string
        full_strings.append(full_string)

    df["alignment_string"] = full_strings

    return df


def generate_match_count_columns(df):

    def count_ones(str, seed=False):
        return str[1:7].count("1") if seed else str.count("1")

    df["pred_num_basepairs"] = df["alignment_string"].apply(count_ones)
    df["pred_seed_basepairs"] = df["alignment_string"].apply(
        count_ones, seed=True)

    return df