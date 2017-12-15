#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""The confdist application."""

import os.path

import csv
import unicodedata
import logging
import argparse
import numpy as np
import pandas as pd


# print log to screen
logging.basicConfig(filename='confdist.log',
                    filemode='w',
                    level=logging.DEBUG, 
                    format='%(levelname)s - %(message)s')
# logging.basicConfig(level=logging.DEBUG, 
#                     format='%(levelname)s - %(message)s')

# read combinations from file
def read_combinations(csv_in):
    """Reads the combinations and their cost from a CSV file and returns a dictionary
    
    The input CSV must have 3 columns in the following order:
    "from", "to", "cost"
    
    All letters are cast to lowercase.
    
    csv_in: string of path to combinations CSV
    """
    with open(csv_in, mode='r', encoding='utf-8') as f_in:
        combinations = {}
        reader = csv.reader(f_in, delimiter=',', quotechar='"')
        next(reader) # skip the header
        for row in reader:
            f = row[0].lower() # from
            to = row[1].lower()
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


def strip_accents_lower(s):
    g = ''.join(c for c in unicodedata.normalize('NFD', s)
                  if unicodedata.category(c) != 'Mn')
    g = g.lower()
    # g.replace() does not alter g itself, it just returns the modified g
    g = g.replace(" ", "")
    return g


def adapted_damerau_levenshtein_distance(string1, string2, combinations):
    """
    Results are rounded up to 2 decimals.
    """
    # for the case of scriptio continua, merge words
    s1 = strip_accents_lower(string1)
    s2 = strip_accents_lower(string2)
    logging.debug("source: %s, target: %s", string1, string2)
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    
    cost_equal = 0
    cost_deletion = 1
    cost_insertion = 1
    cost_substitution = 1
    cost_transposition = 1

    for i in range(-1, lenstr1 + 1):
        d[(i, -1)] = i + 1
    for j in range(-1, lenstr2 + 1):
        d[(-1, j)] = j + 1
    
    for i, li in enumerate(s1):
        for j, lj in enumerate(s2): 
            logging.debug('str1: %s str2: %s', s1[i], s2[j])
            logging.debug("top: %s", d[(i - 1, j)])
            logging.debug("diagonal %s", d[(i - 1, j - 1)])
            logging.debug("left %s", d[(i, j - 1)])
            
            key = "{}-{}".format(s1[i], s2[j])
            
            # for divisions
            if j <= 0:
                key_combo = None
            else:
                combo = s2[j-1] + s2[j]
                key_combo = "{}-{}".format(s1[i], combo)
                logging.debug("division combination: %s", key_combo)
            
            if j <= 1:
                key_prevcombo = None
            else:
                combo = s2[j-2] + s2[j-1]
                key_prevcombo = "{}-{}".format(s1[i], combo)
                logging.debug("division prev combination: %s", key_prevcombo)
            
            # for contractions
            if i <= 0:
                key_c_combo = None
            else:
                combo = s1[i-1] + s1[i]
                key_c_combo = "{}-{}".format(combo, s2[j])
                logging.debug("contraction combination: %s", key_c_combo)
                
            # for complex substituions
            if i <= 0 or j <= 0:
                key_complex = None
            else:
                combo_j = s2[j-1] + s2[j]
                combo_i = s1[i-1] + s1[i]
                key_complex = "{}-{}".format(combo_i, combo_j)
                logging.debug("complex combination: %s", key_complex)
            
            if i > 0 and j > 0:
                logging.debug("i > 0 and j > 0 in source and target")
                if key_complex in combinations.keys():
                    logging.debug("valid complex combination")
                    # cost = cost until the combination + cost of combination
                    d[(i, j)] = d[(i - 2, j - 2)] + combinations[key_complex] 
                elif s1[i] == s2[j]:
                    logging.debug("letters are equal...")
                    if s1[i-1] == s2[j-1]:
                        logging.debug("and also previous are equal -> do not form a combination, take diagonal")
                        d[(i, j)] = d[(i - 1, j - 1)]
                    elif key_combo in combinations.keys():
                        logging.debug("and do form a combination, take diagonal")
                        # take the min. from the previous step, because the total
                        # cost in case of a combination is, the cost of operation
                        # on j-1 plus the cost of the combination
                        # In case of a combination, j counts as if it was j-1,
                        # but including the combination operation
                        d[(i, j)] = min(
                                        d[(i-1, j-1)],
                                        d[(i-1, j-2)],
                                        d[(i, j-2)]
                                        ) + combinations[key_combo]
                    else:
                        logging.debug("previous are unequal and not combination")
                        d[(i, j)] = d[(i - 1, j - 1)]
                elif key_combo in combinations.keys():
                    logging.debug("division combination is valid")
                    # take the min. from the previous step, because the total
                    # cost in case of a combination is, the cost of operation
                    # on j-1 plus the cost of the combination
                    # In case of a combination, j counts as if it was j-1,
                    # but including the combination operation
                    d[(i, j)] = min(
                                    d[(i-1, j-1)],
                                    d[(i-1, j-2)],
                                    d[(i, j-2)]
                                    ) + combinations[key_combo]
                elif key_c_combo in combinations.keys():
                    # the cost is computed as for the previous letter in s1,
                    # but with accounting for the combination
                    d[(i, j)] = min(
                                    d[(i-2, j)],
                                    d[(i-2, j-1)],
                                    d[(i-1, j-1)]
                                    ) + combinations[key_c_combo]
                elif key in combinations.keys():
                    logging.debug("single letter combination")
                    d[(i, j)] = min(
                                    d[(i-1, j-1)],
                                    d[(i-1, j-2)],
                                    d[(i, j-2)]
                                    ) + combinations[key]
                else:
                    logging.debug("not equal and not combination")
                    d[(i, j)] = min(
                           d[(i - 1, j)] + cost_deletion,  # deletion
                           d[(i, j - 1)] + cost_insertion,  # insertion
                           d[(i - 1, j - 1)] + cost_substitution,  # substitution
                           )
            else:
                logging.debug("the first letters in source or target")
                
                if i == 0 and j == 0:
                    logging.debug("fist letters in both")
                    if s1[i] == s2[j]:
                        logging.debug("if equal, take diagonal")
                        d[(i, j)] = d[(i - 1, j - 1)]
                    elif key in combinations.keys():
                        logging.debug("if combination, diagonal plus combination cost")
                        d[(i, j)] = d[(i - 1, j - 1)] + combinations[key]
                    else:
                        # substitution
                        d[(i, j)] = d[(i - 1, j - 1)] + cost_substitution
                elif i == 0 or j == 0:
                    if key_prevcombo in combinations.keys():
                        logging.debug("division prev combination is valid")
                        # then the current cost is the cost until j-1 plus
                        # inserting a new letter
                        d[(i, j)] = d[(i, j-1)] + cost_insertion
                    elif key_combo in combinations.keys():
                        logging.debug("division combination is valid")
                        # the cost is computed as for the previous letter in s2,
                        # but with accounting for the combination
                        d[(i, j)] = min(
                                        d[(i-1, j-1)],
                                        d[(i-1, j-2)],
                                        d[(i, j-2)]
                                        ) + combinations[key_combo]
                    elif key_c_combo in combinations.keys():
                        # the cost is computed as for the previous letter in s1,
                        # but with accounting for the combination
                        d[(i, j)] = min(
                                        d[(i-2, j)],
                                        d[(i-2, j-1)],
                                        d[(i-1, j-1)]
                                        ) + combinations[key_c_combo]
                    else:
                        # not a combination
                        logging.debug("not equal and not combination")
                        d[(i, j)] = min(
                               d[(i - 1, j)] + cost_deletion,  # deletion
                               d[(i, j - 1)] + cost_insertion,  # insertion
                               d[(i - 1, j - 1)] + cost_substitution,  # substitution
                               )
                else:
                    logging.error("unexpected path: i != 0 and j != 0")
                    pass
    distance = round(d[lenstr1 - 1, lenstr2 - 1], 2)
    logging.debug("distance is %s \n", distance)
    return distance


def compare_words(word_list, combinations_dict, verbose=True):
    """Takes a list of dictionaries and compares the words.
    
    The each dictionary must have at least the following keys:
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
        pair["d_computed"] = adapted_damerau_levenshtein_distance(source, target, combinations_dict)
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


def distance_matrix(word_list, combinations_dict):
    """Compares a list of words to each other. Returns a pandas data frame.
    
    https://stackoverflow.com/questions/37428973/string-distance-matrix-in-python
    https://stackoverflow.com/questions/11106536/adding-row-column-headers-to-numpy-matrices
    """
    List1 = word_list
    List2 = word_list
    Matrix = np.zeros((len(List1), len(List2)), dtype=np.float)
    
    for i in range(0, len(List1)):
        for j in range(0, len(List2)):
            Matrix[i, j] = adapted_damerau_levenshtein_distance(List1[i],
                                                             List2[j],
                                                             combinations_dict)
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
        "-c",
        help="The letter confusion table (CSV format).")
    parser.add_argument(
        "-i",
        help="Input table (CSV format).")
    parser.add_argument(
        "-o",
        help="Output table (CSV format).")
    parser.add_argument(
        "-l",
        help="Compare all words in the source column of input table to each other? (yes/no)",
        default="no",
        type=str)

    args = parser.parse_args()
    args_in = {}
    args_in['conf_table'] = os.path.abspath(args.c)
    args_in['file_in'] = os.path.abspath(args.i)
    args_in['file_out'] = os.path.abspath(args.o)
    if "y" in args.l.lower():
        args_in['as_list'] = True
    else:
        args_in['as_list'] = False

    word_list = read_words(args_in['file_in'], as_list=args_in['as_list'])
    conf_table = read_combinations(args_in['conf_table'])
    
    if args_in['as_list']:
        print("Computing distances...")
        m = distance_matrix(word_list, conf_table)
        df = pd.DataFrame(m, index=word_list, columns=word_list)
        print("Writing file...")
        df.to_csv(args_in['file_out'], index=True, header=True, sep=',', 
                  encoding='utf-8')
    else:
        print("Computing distances...")
        a = compare_words(word_list, conf_table, verbose=False)
        print("Writing file...")
        write_results(a, args_in['file_out'])


if __name__ == '__main__':
    main()


