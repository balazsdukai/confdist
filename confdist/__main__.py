#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""The confdist application."""

import os.path

import csv
import unicodedata
import argparse
import numpy as np
import pandas as pd



def read_combinations(csv_in, case):
    """Reads the combinations and their cost from a CSV file and returns a dictionary
    
    The input CSV must have 3 columns in the following order:
    "from", "to", "cost"
    
    Case of all letters can be casted to lowercase, UPPERCASE or can be left as oRIgiNAl.
        
        csv_in: string of path to combinations CSV
    """
    with open(csv_in, mode='r', encoding='utf-8') as f_in:
        combinations = {}
        reader = csv.reader(f_in, delimiter=',', quotechar='"')
        next(reader) # skip the header
        for row in reader:
            if case == "lower":
                f = row[0].lower() # from
                to = row[1].lower()
            elif case == "upper":
                f = row[0].upper() # from
                to = row[1].upper()
            else:
                f = row[0]
                to = row[1]
            cost = float(row[2].replace(',','.')) # replace decimal mark from , to . but this wont work in case of 1,200.30
            key = '-'.join([f, to])
            combinations[key] = cost
            
    return combinations


def read_words(csv_in, as_list=False):
    """Reads word pairs into a dictionary or list suitable for compare_words()
    
    The CSV must have at least "source" and "target" columns if as_list is
    False. If present, the "distance" column is read from the CSV. The
    "distance" column hold the true distances, which are compared by
    compare_words() to the computed distances for validation.
    If as_list is True, only the "source" column is read into a list.
    """
    with open(csv_in, mode='r', encoding='utf-8') as f_in:
        out = []
        reader = csv.DictReader(f_in, delimiter=',', quotechar='"')
        if as_list:
            for row in reader:
                out.append(row["source"])
        else:
            for row in reader:
                try:
                    d = {"source": row["source"],
                         "target": row["target"],
                         "distance": float(row["distance"])}
                except KeyError:
                    d = {"source": row["source"],
                         "target": row["target"]}
                out.append(d)
        return out


def strip_accents_lower(s, case):
    """Strip all diacritics from strings and change the case of the letters to lowercase, UPPERCASE (optional).
        
    Stripping diacritics is important for best emulation of Greek Maiuscule text.
    Case changing is of importance to best fit the provided confusion table.
    """
    g = ''.join(c for c in unicodedata.normalize('NFD', s)
                 if unicodedata.category(c) != 'Mn')
    if case == "lower":
        g = g.lower()
    elif case == "upper":
        g = g.upper()
    else:
        pass
    g = g.replace(" ", "")
    return g


def confusion_distance(string1,
                       string2,
                       combinations,
                       rounding,
                       case,
                       cost_substitution,
                       cost_deletion,
                       cost_insertion
                      ):

    """Compute the confusion distance between a pair of words which can occur in transcription of texts.
        
    This algoritm is an adaptation of an existing implementation of
    the original algoritm by Levensthein.
    
    The existing implementation: 
   
    Parameters
    ----------
    string1 : str
        Source string
    string2 : str
        Target string
    combinations : dict
        Confusion table and costs
    rounding : int
        Number of decimal points for rounding of computed distance (default=2)
    case : string
        Convert source and target to lower/upper case?
    cost_substitution : float
        Cost of substitution (default=1.0)
    cost_deletion : float
        Cost of deletion (default=1.0)
    cost_insertion : float
        Cost of insertion (default=1.0)
    """
    
    # for the case of scriptio continua, merge words
    s1 = strip_accents_lower(string1, case)
    s2 = strip_accents_lower(string2, case)
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)

    """ Create two virtual colums and two virtual rows
        to enable the calculation of contractions, divisions
        and explosions for first characters in either string. 
    """

    for i in range(-2, lenstr1 + 2):
        d[(i, -2)] = i + 1
    for j in range(-2, lenstr2 + 2):
        d[(-2, j)] = j + 1
    
    for i in range(-1, lenstr1 + 1):
        d[(i, -1)] = i + 1
    for j in range(-1, lenstr2 + 1):
        d[(-1, j)] = j + 1
        
    for i, li in enumerate(s1):      
        for j, lj in enumerate(s2):

            """ Set parameters to default values. These values will be used 
                if the compared characters (see below) are not in the confusion 
                table. Otherwise, they will be overwritten by the values for
                the particular combinations in the confusion table. 
            """

            cost_equal = 0
            cost_deletion = 1
            cost_insertion = 1
            cost_substitution = 1
            cost_explosion = 3
            cost_contraction = 3
            cost_complex_substitution = 5
            cost_transposition = 5
            
            """ If present, replace the default parameters
                by the values in the confusion table. 
            """
    
            # for simple | i.e. a-b comparison
            key = "{}-{}".format(s1[i], s2[j])
            cost_substitution = combinations.get(key,cost_substitution)
            
            # for transpositions | i.e a = b-1 and a-1 = b
            if i and j and s1[i]==s2[j-1] and s1[i-1] == s2[j]:
                cost_transposition = 1

            # for explosions | i.e. a-bc comparison
            if j <= 0:
                explosion = None
            else:
                combo = s2[j-1] + s2[j]
                explosion = "{}-{}".format(s1[i], combo)
                cost_explosion = combinations.get(explosion,cost_explosion)
             
            # for contractions | i.e. ab-c comparision
            if i <= 0:
                contraction = None
            else:
                combo = s1[i-1] + s1[i]
                contraction = "{}-{}".format(combo, s2[j])
                cost_contraction = combinations.get(contraction,cost_contraction)

            # for complex substituions | i.e. ab-cd comparision
            if i <= 0 or j <= 0:
                complex_substitution = None
            else:
                combo_j = s2[j-1] + s2[j]
                combo_i = s1[i-1] + s1[i]
                complex_substitution = "{}-{}".format(combo_i, combo_j)
                cost_complex_substitution = combinations.get(complex_substitution,cost_complex_substitution)

            ##################
            ### levensthein
            ##################

            if i == 0 and j == 0:
                """ Both compared characters are on the first 
                    position within the equated words. 
                """

                if s1[i] == s2[j]:
                    d[(i, j)] = d[(i - 1, j - 1)]
                else:
                    d[(i, j)] = d[(i - 1, j - 1)] + cost_substitution            
            
            elif i > 0 or j > 0:
                """ The compared characters are not on the first 
                    position of the equated words. 
                """

                if s1[i] == s2[j]:
                    """ If characters in source and target are equal,
                        the diagonal position should be taken. 
                    """

                    d[(i, j)] = d[(i-1, j-1)]
                    
                else:
                    """ Characters in source and target are dissimmilar,
                        the minimum count of (i-1,j) + the cost of deletion, 
                                             (i-1,j-1) + the cost of substitution,
                                             (i,j-1) + the cost of insertion,
                                             (i-2,j-2)+ the cost of complex substitution,
                                             (i-2,j-1) + the cost of contraction,
                                             (i-1,j-2) + the cost of explosion
                        should be taken. 
                    """

                    d[(i, j)] = min(
                                    d[(i, j-1)] + cost_insertion,            
                                    d[(i-1, j-1)] + cost_substitution,
                                    d[(i-1, j)] + cost_deletion,
                                    d[(i-2,j-2)] + cost_transposition,
                                    d[(i-1, j-2)] + cost_explosion,
                                    d[(i-2, j-1)] + cost_contraction,
                                    d[(i-2, j-2)] + cost_complex_substitution
                                    )
            else:
                pass

    lenstr1 = int(lenstr1)
    lenstr2 = int(lenstr2)
    rounding = int(rounding)
    
    distance = round(d[lenstr1 - 1, lenstr2 - 1], rounding)
    
    return distance


def compare_words(word_list, combinations_dict,
                    verbose,
                    rounding,
                    case,
                    cost_substitution,
                    cost_deletion,
                    cost_insertion):
    """Takes a list of dictionaries and compares the words.
    
    Each dictionary must have at least the following keys:
        "source"
        "target"
    If an expected distance "distance" is provided, it is compared to the
    computed distance and the result printed on the screen.
    Returns word_list with the computed distanced "d_computed" appended to each
    element.
    """
    
    sample = word_list
    for pair in sample:
        source = pair["source"]
        target = pair["target"]
        pair["d_computed"] = confusion_distance(source,
                                                target,
                                                combinations_dict,
                                                rounding=rounding,
                                                case=case,
                                                cost_substitution=cost_substitution,
                                                cost_deletion=cost_deletion,
                                                cost_insertion=cost_insertion)
    if verbose:
        try:
            diff = abs(pair["d_computed"] - pair["distance"]) < 0.0001
            msg = "{}; the distance between {} and {} is {}".format(diff, source,
                                                            target, pair["d_computed"])
            print(msg)
        except KeyError:
            pass
    else:
        pass
    
    return sample


def distance_matrix(word_list, combinations_dict,
                    rounding,
                    case,
                    cost_substitution,
                    cost_deletion,
                    cost_insertion):
    """Compares a list of words to each other. Returns a pandas data frame.
    
    https://stackoverflow.com/questions/37428973/string-distance-matrix-in-python
    https://stackoverflow.com/questions/11106536/adding-row-column-headers-to-numpy-matrices
    """
    
    List1 = word_list
    List2 = word_list
    Matrix = np.zeros((len(List1), len(List2)), dtype=np.float)

    for i in range(0, len(List1)):
        for j in range(0, len(List2)):
            Matrix[i, j] = confusion_distance(List1[i],
                                              List2[j],
                                              combinations_dict,
                                              rounding=rounding,
                                              case=case,
                                              cost_substitution=cost_substitution,
                                              cost_deletion=cost_deletion,
                                              cost_insertion=cost_insertion)

    return Matrix


def write_results(res, file_out):
    """Write a dictionary of word pairs and distances into a CSV.
    
    For the distance matrix, use the pandas to_csv() function.
    """
    with open(file_out, mode='w', encoding='utf-8') as f_out:
        fieldnames = ["source", "target", "distance", "d_computed"]
        writer = csv.DictWriter(f_out, fieldnames=fieldnames,
                                 delimiter=',', quotechar='"')
        writer.writeheader()
        for rowdict in res:  
            writer.writerow(rowdict)


def main():
    parser = argparse.ArgumentParser(description="Compute transcriptional edit distance between pairs of words.")
    parser.add_argument(
        "confusion_table",
        help="The letter confusion table (CSV format).")
    parser.add_argument(
        "input_file",
        help="Input table (CSV format).")
    parser.add_argument(
        "output_file",
        help="Output table (CSV format).")
    parser.add_argument(
        "-l",
        help="Compare all words in the source column of input table to each other? (yes/no) (default=no)",
        default="no",
        type=str)
    parser.add_argument(
        "-rounding",
        help="Number of decimal points for rounding of computed distance (default=2)",
        default=2,
        type=int)
    parser.add_argument(
        "-case",
        help="Convert source and target to lower/upper case, or leave as original? (lower/upper/original) (default=lower)",
        default="lower",
        type=str)
    parser.add_argument(
        "-cs",
        help="Cost of substitution (default=1.0)",
        default=1.0,
        type=float)
    parser.add_argument(
        "-cd",
        help="Cost of deletion (default=1.0)",
        default=1.0,
        type=float)
    parser.add_argument(
        "-ci",
        help="Cost of insertion (default=1.0)",
        default=1.0,
        type=float)
    parser.add_argument(
        "-v",
        help="Verbose mode (yes/no) (default=yes)",
        default="yes",
        type=str)
    
    args = parser.parse_args()
    args_in = {}
    args_in['conf_table'] = os.path.abspath(args.confusion_table)
    args_in['file_in'] = os.path.abspath(args.input_file)
    args_in['file_out'] = os.path.abspath(args.output_file)
    args_in['rounding'] = args.rounding
    args_in['cs'] = args.cs
    args_in['cd'] = args.cd
    args_in['ci'] = args.ci
    args_in['case'] = args.case.lower()
    if "y" in args.l.lower():
        args_in['as_list'] = True
    else:
        args_in['as_list'] = False
    if "y" in args.v.lower():
        args_in['verbose'] = True
    else:
        args_in['verbose'] = False
        
        
    word_list = read_words(args_in['file_in'], as_list=args_in['as_list'])
    conf_table = read_combinations(args_in['conf_table'], args_in['case'])
    
    if args_in['as_list']:
        print("Computing distances...")
        m = distance_matrix(word_list, conf_table,
                            rounding=args_in['rounding'],
                            case=args_in['case'],
                            cost_deletion=args_in['cd'],
                            cost_insertion=args_in['ci'],
                            cost_substitution=args_in['cs'])
        df = pd.DataFrame(m, index=word_list, columns=word_list)
        print("Writing file...")
        df.to_csv(args_in['file_out'], index=True, header=True, sep=',',
                          encoding='utf-8')
    else:
        print("Computing distances...")
        a = compare_words(word_list, conf_table,
                            verbose=args_in['verbose'],
                            rounding=args_in['rounding'],
                            case=args_in['case'],
                            cost_deletion=args_in['cd'],
                            cost_insertion=args_in['ci'],
                            cost_substitution=args_in['cs'])
    print("Writing file...")
    write_results(a, args_in['file_out'])


