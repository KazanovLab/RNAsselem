import sys
from math import sqrt
from operator import itemgetter
from os.path import isfile
from os.path import basename


index_column_number = 0
connection_column_number = 4
symbol_column_number = 6


# public
def load_wuss(wuss_path):

    if isinstance(wuss_path, list):
        return wuss_path

    else:
        input_file = open(wuss_path, "r")
        input_file.readline()
        genome = []

        for line in input_file:
            nucleotide = line.split()

            for i in range(len(nucleotide)):
                symbol = nucleotide[i]
                if symbol.isdecimal():
                    nucleotide[i] = int(symbol)

            genome.append(nucleotide)
        return genome


def index_minus_one(wuss):

    for i in range(len(wuss)):
        nucleotide = wuss[i]

        if isinstance(nucleotide, list):
            for j in range(len(nucleotide)):
                element = nucleotide[j]
                if isinstance(element, int):
                    wuss[i][j] = element - 1

        else:
            if isinstance(nucleotide, int):
                wuss[i] = nucleotide - 1

    return wuss


# private
def get_columns(wuss_list, column_numbers):

    if isinstance(column_numbers, int):
        column_numbers = [column_numbers]
    new_wuss_list = []
    for nucleotide in wuss_list:
        only_necessary = []
        for i in column_numbers:
            only_necessary.append(nucleotide[i])
        if len(only_necessary) == 1:
            new_wuss_list.append(only_necessary[0])
        else:
            new_wuss_list.append(only_necessary)

    return new_wuss_list


# public
def get_total_length(wuss):

    wuss_array = load_wuss(wuss)
    return len(wuss_array)


# public
def get_unpaired(wuss):

    wuss_list = load_wuss(wuss)
    connection_indexes = get_columns(wuss_list, connection_column_number)
    zero = connection_indexes.count(0)
    return zero


# public
def get_paired(wuss):

    total_length = get_total_length(wuss)
    unpaired = get_unpaired(wuss)
    paired = total_length - unpaired
    return paired


# private
def statistics(results_list):

    length_list = []
    for element in results_list:
        length = element["length"]
        length_list.append(length)

    count = len(length_list)
    length_list.sort()
    mean_length = sum(length_list) / count

    if count % 2 == 0:
        r = count // 2
        l = r - 1
        median_length = (length_list[r] + length_list[l]) / 2
    else:
        i = (count + 1) // 2
        median_length = length_list[i]

    sq = []
    for n in length_list:
        m = n - mean_length
        sq.append(m ** 2)
    d = sum(sq) / count
    std_length = sqrt(d)

    total_length = sum(length_list)

    values = [count, round(mean_length, 2), round(median_length, 2), round(std_length, 2), total_length]
    keys = ["count", "mean_length", "median_length", "std_length", "total_length"]
    stat_dict = dict(zip(keys, values))
    return stat_dict


# private
def make_indexes_dict(length, start_index, end_index):

    d = {"length": length, "start_index": start_index + 1, "end_index": end_index + 1}
    return d


# public
def get_external_loops_list(wuss):

    wuss_list = load_wuss(wuss)
    wuss_symbols = get_columns(wuss_list, symbol_column_number)

    external_loops = []
    length = 0

    for i in range(len(wuss_symbols)):
        nucleotide = wuss_symbols[i]
        if nucleotide == ":":
            if length == 0:
                l_index = i + 1
            length += 1
        else:
            if length != 0:
                r_index = i
                external_loops.append({"start_index": l_index, "end_index": r_index, "length": length})
                length = 0

    if length != 0:
        r_index = len(wuss_symbols)
        external_loops.append({"start_index": l_index, "end_index": r_index, "length": length})

    return external_loops


# private
def define_input(input_object):

    if isinstance(input_object, str):
        return "path"
    else:
        first_element = input_object[0]
        if isinstance(first_element, list):
            return "input_list"
        else:
            return "results_list"


# public
def get_external_loops_stat(wuss):

    input_type = define_input(wuss)

    if input_type == "path" or input_type == "input_list":
        wuss = get_external_loops_list(wuss)

    return statistics(wuss)


# public
def get_hairpin_loops_list(wuss):

    wuss_list = load_wuss(wuss)
    wuss_symbols = get_columns(wuss_list, symbol_column_number)

    hairpin_loops = []
    length = 0

    for i in range(len(wuss_symbols)):
        nucleotide = wuss_symbols[i]
        if nucleotide == "_":
            if length == 0:
                l_index = i + 1
            length += 1
        else:
            if length != 0:
                r_index = i
                hairpin_loops.append({"start_index": l_index, "end_index": r_index, "length": length})
                length = 0

    return hairpin_loops


# public
def get_hairpin_loop_stat(wuss):

    input_type = define_input(wuss)

    if input_type == "path" or input_type == "input_list":
        wuss = get_hairpin_loops_list(wuss)

    return statistics(wuss)


# private
def get_bulge_and_internal_loops_list(wuss):

    wuss_list = index_minus_one(load_wuss(wuss))
    genome = get_columns(wuss_list, [index_column_number, connection_column_number, symbol_column_number])
    brackets = ["<", "(", "[", "{"]
    cl_brackets = [">", ")", "]", "}"]

    def interior_or_bulge(left_index, right_index, genome):
        slice = genome[left_index:right_index + 1]
        quantity = 0
        for nucleotide in slice:
            if nucleotide[2] == "-":
                quantity += 1
        return quantity

    length = 0
    interior_loops = []
    bulge_loops = []
    used = []

    for i in range(len(genome)):
        nucleotide = genome[i]
        previous_nucleotide = genome[i - 1] if i != 0 else [[], [], []]

        if nucleotide[2] in brackets:

            if previous_nucleotide[2] == "-":
                right_index = nucleotide[0]
                left_connection_index = genome[right_index][1]
                right_connection_index = genome[left_index][1]
                quantity = interior_or_bulge(left_connection_index, right_connection_index, genome)
                if quantity > 0:
                    total_length = length + quantity
                    dictionary = make_indexes_dict(total_length, left_index + 1, right_connection_index - 1)
                    interior_loops.append(dictionary)
                    used.append([left_connection_index, right_connection_index])
                else:
                    dictionary = make_indexes_dict(length, left_index + 1, right_index - 1)
                    bulge_loops.append(dictionary)
                length = 0

        elif nucleotide[2] in cl_brackets:

            if previous_nucleotide[2] == "-":
                right_index = nucleotide[0]

                if [left_index, right_index] not in used:
                    dictionary = make_indexes_dict(length, left_index + 1, right_index - 1)
                    bulge_loops.append(dictionary)
                length = 0

        elif nucleotide[2] == "-":
            if previous_nucleotide[2] in (brackets + cl_brackets):
                left_index = previous_nucleotide[0]
            length += 1

    return [interior_loops, bulge_loops]


# public
def get_internal_loop_list(wuss):

    result = get_bulge_and_internal_loops_list(wuss)[0]
    return result


# public
def get_internal_loop_stat(wuss):

    input_type = define_input(wuss)

    if input_type == "path" or input_type == "input_list":
        wuss = get_internal_loop_list(wuss)

    return statistics(wuss)


# public
def get_bulge_loop_list(wuss):

    result = get_bulge_and_internal_loops_list(wuss)[1]
    return result


# public
def get_bulge_loop_stat(wuss):

    input_type = define_input(wuss)

    if input_type == "path" or input_type == "input_list":
        wuss = get_bulge_loop_list(wuss)

    return statistics(wuss)


# public
def get_multifurcation_loop_list(wuss):

    wuss_list = load_wuss(wuss)
    genome = index_minus_one(get_columns(wuss_list, [index_column_number, connection_column_number, symbol_column_number]))

    def find_bracket(bracket, array):
        for i in range(len(array)):
            line = array[i]
            if bracket in line:
                return i
                break
        else:
            return -1

    def delete_bracket(array, index):
        beginning = array[index]
        connection_number = beginning[1]
        length = 0
        r = array[index:]

        start_index = -2
        last_element = [-2] * 3

        for nucleotide in r:
            if nucleotide[2] == ",":
                start_index = nucleotide[0]
                break

        for i in range(len(r)):
            nucleotide = r[i]
            if nucleotide[2] == ",":
                length += 1
                last_element = nucleotide
            if nucleotide[0] == connection_number:
                l_index = index + i
                end_index = last_element[0]
                break

        l = array[:index]
        r = array[l_index + 1:]
        new_array = l + r

        multifurcation_loop = make_indexes_dict(length, start_index, end_index)
        return [new_array, multifurcation_loop]

    new_genome = genome
    multifurcation_loops = []

    for bracket in ["<", "(", "[", "{"]:
        while True:
            bracket_index = find_bracket(bracket, new_genome)
            if bracket_index < 0:
                break
            new_genome, multifurcation_loop = delete_bracket(new_genome, bracket_index)
            if multifurcation_loop["length"] != 0:
                multifurcation_loops.append(multifurcation_loop)

    multifurcation_loops.sort(key=itemgetter("start_index"))

    return multifurcation_loops


# public
def get_multifurcation_loops_stat(wuss):

    input_type = define_input(wuss)

    if input_type == "path" or input_type == "input_list":
        wuss = get_multifurcation_loops_list(wuss)

    return statistics(wuss)


# public
def get_stem_list(wuss):

    wuss_list = load_wuss(wuss)
    genome = index_minus_one(get_columns(wuss_list, [index_column_number, connection_column_number, symbol_column_number]))

    brackets = ("(", "[", "<", "{")
    stem_list = []
    length = 0

    for i in range(len(genome)):
        nucleotide = genome[i]

        if nucleotide[2] in brackets:

            previous_nucleotide = genome[i - 1] if i != 0 else [-1] * 3
            if previous_nucleotide[2] not in brackets:
                start_index = i
                end_index = nucleotide[1]
            length += 1
        else:

            if length != 0:
                d = make_indexes_dict(length, start_index, end_index)
                stem_list.append(d)
                length = 0

    return stem_list


# public
def get_stem_stat(wuss):

    input_type = define_input(wuss)

    if input_type == "path" or input_type == "input_list":
        wuss = get_stem_list(wuss)

    stat = statistics(wuss)
    stat["total_length"] *= 2

    return stat


# public
def save_stat(wuss, output_path, is_append=False):

    lines = ["External_loops", "Hairpin_loops", "Internal_loops", "Bulge_loops", "Multifurcation_loops", "Stems"]
    stat = [get_external_loops_stat(wuss), get_hairpin_loops_stat(wuss), get_internal_loop_stat(wuss),
            get_bulge_loops_stat(wuss), get_multifurcation_loops_stat(wuss), get_stem_stat(wuss)]

    if isfile(output_path) and not is_append:
        print("This file is already exist. Delete this file and try again.")
        return -1

    output_file = open(output_path, "a")
    if is_append:
        output_file.write("\n" * 2)

    line1 = basename(wuss)
    output_file.write(line1 + "\n")
    line2 = "\tCount\tMean_length\tMedian_length\tStd_length\tTotal_length\n"
    output_file.write(line2)

    for i in range(6):
        line = lines[i] + "\t"
        stat_values = list(stat[i].values())
        for i in range(5):
            stat_values[i] = str(stat_values[i])
        line += "\t".join(stat_values) + '\n'
        output_file.write(line)
    return 0







