import copy
import sys
from math import sqrt
from operator import itemgetter
from os.path import isfile
from os.path import basename
from copy import deepcopy

index_column_number = 0
letter_column_number = 1
connection_column_number = 4
symbol_column_number = 6
structure_column_number = 7


def ct2ss(ctPath):
    try:
        ctfile = open(ctPath)
    except FileNotFoundError:
        print("File not found.")
    except PermissionError:
        print("Permission error. Cannot open the file.")
    except Exception as e:
        print("An error occurred:", e)

    ct = [0]
    ct2 = [0]
    line = ctfile.readline()
    flds = line.split()
    numlines = int(flds[0])
    i = 1
    line = ctfile.readline()
    while (line):
        flds = line.split()
        if len(flds) != 6:
            print("Error: data line does not contain six columns.")
            sys.exit()
        if flds[0] != str(i):
            print("Error: wrong line number " + flds[0] + " for " + str(i) + "-th line.")
            sys.exit()
        ct.append(int(flds[4]))
        ct2.append(int(flds[4]))
        line = ctfile.readline()
        i += 1

    if numlines != (i - 1):
        print(
            "Error: number of actual lines " + str(i - 1) + " does not equal specified in file " + str(numlines) + ".")
        sys.exit()

    n = numlines

    npairs = 0
    for i in range(1, n + 1):
        if ct[i] != 0 and i < ct[i]:
            npairs += 1

    numpseudo = 0
    rb = [-1] * 26
    ss = [":"] * n

    npairsActual = 0
    stack = []
    stackLoops = []
    stackPseudo = []
    for i in range(1, n + 1):
        if ct[i] == 0:
            stack.append(i)
        elif ct[i] > i:
            stack.append(i)
        else:
            ispair = 0
            numstr = 0
            minlevel = -1
            while len(stack) > 0:
                j = stack.pop()
                if j < 0:
                    numstr += 1
                    if j < minlevel:
                        minlevel = j
                elif ct[j] == i:
                    ispair = 1
                    npairsActual += 1
                    if numstr > 1 and minlevel > -4:
                        minlevel -= 1
                    if minlevel == -1:
                        ss[j - 1] = '<'
                        ss[i - 1] = '>'
                    elif minlevel == -2:
                        ss[j - 1] = '('
                        ss[i - 1] = ')'
                    elif minlevel == -3:
                        ss[j - 1] = '['
                        ss[i - 1] = ']'
                    elif minlevel == -4:
                        ss[j - 1] = '{'
                        ss[i - 1] = '}'
                    stack.append(minlevel)
                    # loops
                    while len(stackLoops) > 0:
                        j = stackLoops.pop()
                        if numstr == 0:
                            ss[j - 1] = '_'
                        elif numstr == 1:
                            ss[j - 1] = '-'
                        elif numstr > 1:
                            ss[j - 1] = ','
                    break
                elif ct[j] == 0:
                    if ct2[j] == 0:
                        stackLoops.append(j)
                else:
                    stackPseudo.append(j)

            if ispair == 0:
                print("Cannot find left partner " + ct[i] + " of base " + i)
                sys.exit()

        if len(stackPseudo) > 0:
            k = 0
            lbound = ct[i]
            rbound = lbound + 1
            xpk = -1

            while len(stackPseudo) > 0:
                j = stackPseudo.pop()
                print(k)
                k = rbound - 1
                while k > lbound:
                    if ct[k] == 0:
                        k -= 1
                        continue
                    elif ct[k] > rbound:
                        k -= 1
                        continue
                    elif ct[k] == j:
                        break
                    else:
                        k = lbound
                        break
                    k -= 1

                print(k)
                print(lbound)
                print(rbound)
                if k == lbound:
                    numpseudo += 1
                    xpk += 1
                    while j < rb[xpk]:
                        xpk += 1
                    if rbound < ct[j]:
                        lbound = rbound
                    else:
                        lbound = ct[i]
                    rbound = ct[j]

                npairsActual += 1

                if (xpk + ord('a')) <= ord('z'):
                    if ct[j] > rb[xpk]:
                        rb[xpk] = ct[j]
                    print(xpk)
                    ss[j - 1] = chr(xpk + ord('A'))
                    ss[ct[j] - 1] = chr(xpk + ord('a'))

                    ct[j] = 0
                    ct[ct2[j]] = 0

                else:
                    print("Too many pseudoknots to describe by letters.")
                    sys.exit()

    return (ss)


def ct2ctwuss(ctPath, outPath):
    """ Convert connectivity table (CT) format to CT-modified format with WUSS notation.
    
    Function converts connectivity table (CT) file with RNA secondary structure into the same format with WUSS annotation in additional column.
    
    Parameters:
        ctPath (str): path to the connectivity table (CT) file
        outPath (ctr): path for creating output CT-WUSS file
    """

    ss = ct2ss(ctPath)

    fr = open(ctPath)
    fo = open(outPath, "w")
    # skip header
    header = fr.readline()
    fo.write(header)
    s = fr.readline()
    i = 0
    while (s):
        s = s.rstrip() + ' ' + ss[i] + '\n'
        fo.write(s)
        s = fr.readline()
        i += 1
    fr.close()
    fo.close()


def ct2wuss(ctPath, outPath):
    """ Convert connectivity table (CT) format to WUSS format.
        
        Function converts connectivity table (CT) file with RNA secondary structure into the WUSS annotation file.
        
        Parameters:
            ctPath (str): path to the connectivity table (CT) file
            outPath (ctr): path for creating output WUSS file
        """

    ss = ct2ss(ctPath)

    fr = open(ctPath)
    fo = open(outPath, "w")
    # skip header
    head = fr.readline()
    flds = head.split()
    headName = flds[len(flds) - 1]
    s = fr.readline()
    i = 0
    seq = ''
    while (s):
        seq += s.split()[1]
        s = fr.readline()

    fo.write(">" + headName + '\n')
    fo.write(seq + '\n')
    fo.write(''.join(ss) + '\n')

    fr.close()
    fo.close()


# public
def get_sequences(ctwuss, intervals, positon=-2):
    """ Returns the list of nucleotide sequences for the given list of intervals.

        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
            intervals (list): list of pairs reflecting start/end positions of given intervals

        Returns:
            sequences (list): list of nucleotide sequences
    """

    wuss_list = load_ctwuss(ctwuss)
    sequences = []
    for interval in intervals:
        sequence = ""
        l = interval[0] - 1
        r = interval[1] - 1

        for i in range(l, r + 1):
            letter = wuss_list[i][letter_column_number]
            if i == positon - 1:
                sequence += letter
            else:
                sequence += letter.lower()

        sequences.append(sequence)

    return sequences


#private
def get_structure_dict(wuss_list, length, intervals):
    sequences = get_sequences(wuss_list, intervals)
    structure_dict = {"length": length, "intervals": intervals, "sequences": sequences}
    return structure_dict


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


# private
def change_letters(wuss_list, structure_column_index):
    wuss_list_copy = copy.deepcopy(wuss_list)

    def change_structure(index, wuss_list_copy):
        structures = ["hairpin_loop", "external_loop", "multifurcation_loop",
                      "indefinite_loop", "internal_loop", "bulge_loop"]
        previous_nucleotide = wuss_list_copy[index - 1][structure_column_index]
        next_nucleotide = wuss_list_copy[index + 1][structure_column_index]

        if previous_nucleotide in structures:
            appropriate_structure = previous_nucleotide
        elif next_nucleotide in structures:
            appropriate_structure = next_nucleotide
        else:
            appropriate_structure = change_structure(index + 1, wuss_list_copy)

        wuss_list_copy[index][structure_column_index] = appropriate_structure
        return appropriate_structure

    for i in range(len(wuss_list_copy)):
        if wuss_list_copy[i][structure_column_index] == "pseudoknot":
            change_structure(i, wuss_list_copy)

    return wuss_list_copy


# private
def get_bulge_and_internal_loop_list(wuss_list):

    brackets = ["<", "(", "[", "{"]
    cl_brackets = [">", ")", "]", "}"]
    loops = ["bulge_loop", "internal_loop", "indefinite_loop"]
    length = 0
    interior_loops = []
    bulge_loops = []
    used = []
    left_index = -1

    for i in range(len(wuss_list)):
        nucleotide = wuss_list[i]
        previous_nucleotide = wuss_list[i - 1] if i != 0 else 8*[-3]

        if nucleotide[symbol_column_number] in brackets:

            if previous_nucleotide[structure_column_number] in loops:
                right_index = nucleotide[index_column_number] - 1
                left_connection_index = wuss_list[right_index][connection_column_number] - 1
                right_connection_index = wuss_list[left_index][connection_column_number] - 1
                quantity = right_connection_index - left_connection_index - 1
                if quantity > 0:
                    total_length = length + quantity
                    intervals = [[left_index + 2, right_index], [left_connection_index + 2, right_connection_index]]
                    d = get_structure_dict(wuss_list, total_length, intervals)
                    interior_loops.append(d)
                    used.append([left_connection_index, right_connection_index])
                else:
                    dictionary = get_structure_dict(wuss_list, length, [[left_index + 2, right_index]])
                    bulge_loops.append(dictionary)
                length = 0

        elif nucleotide[symbol_column_number] in cl_brackets:

            if previous_nucleotide[structure_column_number] in loops:
                right_index = nucleotide[index_column_number] - 1

                if [left_index, right_index] not in used:
                    dictionary = get_structure_dict(wuss_list, length, [[left_index + 2, right_index]])
                    bulge_loops.append(dictionary)
                length = 0

        elif nucleotide[structure_column_number] in loops:
            if previous_nucleotide[symbol_column_number] in (brackets + cl_brackets):
                left_index = previous_nucleotide[index_column_number] - 1
            length += 1

    return [interior_loops, bulge_loops]


# public
def load_ctwuss(ctwuss_path):
    """ Loads CT-modified file with WUSS annotation in additional column.
        
        Parameters:
            ctwuss_path (str): path to the connectivity table (CT) file
        
        Returns:
            genome (list[list]): list of nucleotides with pairing info
    """

    if isinstance(ctwuss_path, list):
        return ctwuss_path

    else:
        input_file = open(ctwuss_path, "r")
        input_file.readline()
        wuss_list = []

        for line in input_file:
            nucleotide = line.split()

            for i in range(len(nucleotide)):
                symbol = nucleotide[i]
                if symbol.isdecimal():
                    nucleotide[i] = int(symbol)

            symbol = nucleotide[-1]

            if symbol == ":":
                structure_type = "external_loop"
            elif symbol == ",":
                structure_type = "multifurcation_loop"
            elif symbol == "_":
                structure_type = "hairpin_loop"
            elif symbol == "-":
                structure_type = "indefinite_loop"
            elif symbol in ["<", "(", "[", "{", ">", ")", "]", "}"]:
                structure_type = "stem"
            else:
                structure_type = "pseudoknot"

            nucleotide.append(structure_type)
            wuss_list.append(nucleotide)

        wuss_list = change_letters(wuss_list, symbol_column_number + 1)
        bulge_loops = get_bulge_and_internal_loop_list(wuss_list)[1]
        bulge_indexes = []

        for bulge_loop in bulge_loops:
            interval = bulge_loop["intervals"][0]
            interval_range = list(range(interval[0], interval[1] + 1))
            bulge_indexes += interval_range

        for i in range(len(wuss_list)):
            nucleotide = wuss_list[i]
            if nucleotide[structure_column_number] == "indefinite_loop":
                if i + 1 in bulge_indexes:
                    wuss_list[i][structure_column_number] = "bulge_loop"
                else:
                    wuss_list[i][structure_column_number] = "internal_loop"

        return wuss_list

#private
def pseudoknot_to_symbol(ctwuss_path, output_name):
    wuss_list = load_ctwuss(ctwuss_path)
    output_file = open(output_name, "w")
    input_file = open(ctwuss_path, "r")
    name = input_file.readline()
    output_file.write(name)
    for i in wuss_list:
        line = input_file.readline()
        s = i[structure_column_number]
        new = i[symbol_column_number]

        if s == "external_loop":
            new = ":"
        elif s == "multifurcation_loop":
            new = ","
        elif s == "hairpin_loop":
            new = "_"
        elif s == "bulge_loop" or s == "internal_loop":
            new = "-"

        line = line[:-2] + new + "\n"
        output_file.write(line)

    return 0


# private
def get_interval(structures_list, position_index):

    break_flag = False
    for structure_dictionary in structures_list:
        intervals = structure_dictionary["intervals"]

        for interval in intervals:
            if interval[0] <= position_index <= interval[1]:
                break_flag = True
                break

        if break_flag:
            break

    return structure_dictionary


# public
def get_total_length(wuss):
    wuss_array = load_ctwuss(wuss)
    return len(wuss_array)


# public
def get_unpaired(ctwuss):
    """ Returns the number of unpaired nucleotides in genome.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            unpaired_cnt (int): number of unpaired nucleotides
    """

    wuss_list = load_ctwuss(ctwuss)
    connection_indexes = get_columns(wuss_list, connection_column_number)
    unpaired_cnt = connection_indexes.count(0)
    return unpaired_cnt


# public
def get_paired(ctwuss):
    """ Returns the number of paired nucleotides in genome.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            paired_cnt (int): number of paired nucleotides
    """

    total_length = get_total_length(ctwuss)
    unpaired = get_unpaired(ctwuss)
    paired_cnt = total_length - unpaired
    return paired_cnt


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


# public
def get_external_loop_list(ctwuss):
    """ Returns a list of external loops with a dictionary of their properties.
               
       Parameters:
        ctwuss (str/list): path to the CT-modified file or loaded file

       Returns:
            external_loops (list[dict]): list of the external loop properties, including:
                - length - element length,
                - intervals - list of the start/stop positions,
                - sequences - DNA sequences at the intervals.
    """

    wuss_list = load_ctwuss(ctwuss)
    external_loops = []                                                                      
    length = 0

    for i in range(len(wuss_list)):
        nucleotide = wuss_list[i][structure_column_number]

        if nucleotide == "external_loop":
            if length == 0:
                l_index = i + 1
            length += 1
        else:
            if length != 0:
                r_index = i
                d = get_structure_dict(wuss_list, length, [[l_index, r_index]])
                external_loops.append(d)
                length = 0

    if length != 0:
        r_index = len(wuss_list)
        d = get_structure_dict(wuss_list, length, [[l_index, r_index]])
        external_loops.append(d)

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
def get_external_loop_stat(ctwuss):
    """ Returns statistics on external loops in genome.
        
        Resulting dictionary contains following statistical keys:
        count: number of structural elements
        mean_length: average length of the structural element
        median_length: median length of the structural element
        std_length: standard deviation of the length of the structural elements
        total_length: total length of the structural elements in genome
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file

        Returns:
            stat (dict): dictionary with statistics on external loops
    """

    input_type = define_input(ctwuss)

    if input_type == "path" or input_type == "input_list":
        elements = get_external_loop_list(ctwuss)
    else:
        elements = deepcopy(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_stem_list(ctwuss):
    """ Returns a list of stems with a dictionary of their properties.

        Parameters:
            ctwuss (str): path to the CT-modified file or loaded file

        Returns:
            stem_list (list[dict]): list of the with stem properties, including:
             - length - element length,
             - intervals - list of the start/stop positions,
             - sequences - DNA sequences at the intervals.

    """

    wuss_list = load_ctwuss(ctwuss)
    pairs = [["(", ")"], ["[", "]"], ["<", ">"], ["{", "}"]]
    stem_list = []
    length = 0
    left_index = -1

    intervals = []
    for pair in pairs:
        bracket = pair[0]
        cl_bracket = pair[1]
        for i in range(len(wuss_list)):
            nucleotide = wuss_list[i]
            
            if nucleotide[symbol_column_number] == bracket:
                previous_nucleotide = wuss_list[i - 1]
                next_nucleotide = wuss_list[i + 1] if i != len(wuss_list) - 1 else 8*[-3]
                if length == 0:
                    end_index = nucleotide[connection_column_number] - 1
                length += 1
                if previous_nucleotide[symbol_column_number] != bracket:
                    left_index = nucleotide[index_column_number] - 1
                if next_nucleotide[symbol_column_number] != bracket:
                    intervals.append([left_index + 1, nucleotide[index_column_number]])

            if nucleotide[symbol_column_number] == cl_bracket:
                previous_nucleotide = wuss_list[i - 1]
                next_nucleotide = wuss_list[i + 1] if i != len(wuss_list) - 1 else 8*[-3]
                length += 1
                if previous_nucleotide[symbol_column_number] != cl_bracket:
                    left_index = nucleotide[index_column_number] - 1
                if next_nucleotide[symbol_column_number] != cl_bracket:
                    intervals.append([left_index + 1, nucleotide[index_column_number]])
                if nucleotide[index_column_number] - 1 == end_index:
                    d = get_structure_dict(wuss_list, length, intervals)
                    stem_list.append(d)
                    length = 0
                    intervals = []

    stem_list.sort(key=lambda x: x["intervals"][0][0])

    return stem_list


# public
def get_hairpin_loop_list(ctwuss):
    """ Returns a list of hairpin loops with a dictionary of their properties.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            hairpin_loops (list[dict]): list of the hairpin loop properties, including:
              - length - element length,
              - intervals - list of the start/stop positions,
              - sequences - DNA sequences at the intervals,
              - stem_length - length of the hairpin's stem.
    """

    wuss_list = load_ctwuss(ctwuss)
    stems = get_stem_list(wuss_list)
    hairpin_loops = []
    length = 0

    for i in range(len(wuss_list)):
        nucleotide = wuss_list[i][structure_column_number]
        if nucleotide == "hairpin_loop":
            if length == 0:
                l_index = i + 1
            length += 1
        else:
            if length != 0:
                r_index = i
                stem_dictionary = get_interval(stems, l_index - 1)
                stem_length = stem_dictionary["length"]
                d = get_structure_dict(wuss_list, length, [[l_index, r_index]])
                d["stem_length"] = stem_length
                hairpin_loops.append(d)
                length = 0

    return hairpin_loops


def get_hairpin_loop_stat(ctwuss):
    """ Returns statistics on hairpin loops in genome.
        
        Resulting dictionary contains following statistical keys:
        count: number of structural elements
        mean_length: average length of the structural element
        median_length: median length of the structural element
        std_length: standard deviation of the length of the structural elements
        total_length: total length of the structural elements in genome
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            stat (dict): dictionary with statistics on hairpin loops
    """

    input_type = define_input(ctwuss)

    if input_type == "path" or input_type == "input_list":
        elements = get_hairpin_loop_list(ctwuss)
    else:
        elements = deepcopy(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_internal_loop_list(ctwuss):
    """ Returns a list of internal loops with a dictionary of their properties.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            interior_loops (list[dict]): list of the internal loop properties, including:
             - length - element length
             - intervals - list of the start/stop positions
             - sequences - DNA sequences at the intervals
    """

    wuss_list = load_ctwuss(ctwuss)
    interior_loops = get_bulge_and_internal_loop_list(wuss_list)[0]
    return interior_loops


# public
def get_internal_loop_stat(ctwuss):
    """ Returns statistics on internal loops in genome.
        
        Resulting dictionary contains following statistical keys:
        count: number of structural elements
        mean_length: average length of the structural element
        median_length: median length of the structural element
        std_length: standard deviation of the length of the structural elements
        total_length: total length of the structural elements in genome
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            stat (dict): dictionary with statistics on internal loops
    """

    input_type = define_input(ctwuss)

    if input_type == "path" or input_type == "input_list":
        elements = get_internal_loop_list(ctwuss)
    else:
        elements = deepcopy(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_bulge_loop_list(ctwuss):
    """ Returns a list of bulge loops with a dictionary of their properties.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            bulges (list[dict]): list of the bulge properties, including:
             - length - element length
             - intervals - list of the start/stop positions
             - sequences - DNA sequences at the intervals
    """

    wuss_list = load_ctwuss(ctwuss)
    bulges = get_bulge_and_internal_loop_list(wuss_list)[1]
    return bulges


# public
def get_bulge_loop_stat(ctwuss):
    """ Returns statistics on bulges in genome.
        
        Resulting dictionary contains following statistical keys:
        count: number of structural elements
        mean_length: average length of the structural element
        median_length: median length of the structural element
        std_length: standard deviation of the length of the structural elements
        total_length: total length of the structural elements in genome
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            stat (dict): dictionary with statistics on bulges
    """

    input_type = define_input(ctwuss)

    if input_type == "path" or input_type == "input_list":
        elements = get_bulge_loop_list(ctwuss)
    else:
        elements = deepcopy(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_multifurcation_loop_list(ctwuss):
    """ Returns a list of multifurcation loops with a dictionary of their properties.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            multifurcation_loops (list[dict]): list of the multifurcation loop properties, including:
             - length - element length,
             - intervals - list of the start/stop positions,
             - sequences - DNA sequences at the intervals.

    """

    wuss_list = load_ctwuss(ctwuss)
    wuss_list_copy = deepcopy(wuss_list)

    def find_bracket(bracket, array):
        for i in range(len(array)):
            line = array[i]
            if bracket in line:
                return i
                break
        else:
            return -1

    def delete_bracket(array, bracket_index):

        beginning = array[bracket_index]
        connection_number = beginning[connection_column_number] - 1
        r = array[bracket_index:]
        start_index = 1
        previous_index = -1

        m_nucleotides = []

        for i in range(len(r)):
            nucleotide = r[i]
            if nucleotide[structure_column_number] == "multifurcation_loop":
                m_nucleotides.append(nucleotide)

            if nucleotide[index_column_number] - 1 == connection_number:
                last_index = bracket_index + i
                break

        if len(m_nucleotides) == 1:
            intervals = [[m_nucleotides[0][0], m_nucleotides[0][0]]]
        else:
            intervals = []
            for i in range(len(m_nucleotides)):

                nucleotide = m_nucleotides[i]
                index = nucleotide[index_column_number] - 1

                if i == 0:
                    start_index = index
                else:
                    if index - previous_index != 1:
                        intervals.append([start_index + 1, previous_index + 1])
                        start_index = index
                    if i == len(m_nucleotides) - 1:
                        intervals.append([start_index + 1, index + 1])

                previous_index = index

        del array[bracket_index: last_index + 1]
        multifurcation_loop = get_structure_dict(wuss_list, len(m_nucleotides), intervals)

        return [array, multifurcation_loop]

    multifurcation_loops = []
    for bracket in ["<", "(", "[", "{"]:
        while True:
            bracket_index = find_bracket(bracket, wuss_list_copy)
            if bracket_index < 0:
                break
            wuss_list_copy, multifurcation_loop = delete_bracket(wuss_list_copy, bracket_index)
            if multifurcation_loop["length"] != 0:
                multifurcation_loops.append(multifurcation_loop)
    multifurcation_loops.sort(key=lambda x: x["intervals"][0][0])

    return multifurcation_loops


# public
def get_multifurcation_loop_stat(ctwuss):
    """ Returns statistics on multifurcation loops in genome.
        
        Resulting dictionary contains following statistical keys:
        count: number of structural elements
        mean_length: average length of the structural element
        median_length: median length of the structural element
        std_length: standard deviation of the length of the structural elements
        total_length: total length of the structural elements in genome
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            stat (dict): dictionary with statistics on multifurcation loops
    """

    input_type = define_input(ctwuss)

    if input_type == "path" or input_type == "input_list":
        elements = get_multifurcation_loop_list(ctwuss)
    else:
        elements = deepcopy(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_stem_stat(ctwuss):
    """ Returns statistics on stems in genome.
        
        Resulting dictionary contains following statistical keys:
        count: number of structural elements
        mean_length: average length of the structural element
        median_length: median length of the structural element
        std_length: standard deviation of the length of the structural elements
        total_length: total length of the structural elements in genome
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
            stat (dict): dictionary with statistics on stems
    """

    input_type = define_input(ctwuss)

    if input_type == "path" or input_type == "input_list":
        elements = get_stem_list(ctwuss)
    else:
        elements = deepcopy(ctwuss)

    return statistics(elements)


#public
def get_pseudoknot_list(ctwuss):

    wuss_list = load_ctwuss(ctwuss)
    wuss_list_copy = deepcopy(wuss_list)
    pseudoknot_dict = {}
    previous_symbol = "0"
    start_index = -1

    for nucleotide in wuss_list_copy:

        symbol = nucleotide[symbol_column_number].lower()

        if not symbol.isalpha():

            if previous_symbol.isalpha():
                interval = [start_index, nucleotide[index_column_number] - 1]
                pseudoknot_dict[previous_symbol]["intervals"].append(interval)

        else:
            if symbol not in pseudoknot_dict:
                pseudoknot_dict[symbol] = {"length": 0, "intervals": []}

            if previous_symbol.isalpha():

                if previous_symbol != symbol:
                    interval = [start_index, nucleotide[index_column_number] - 1]
                    pseudoknot_dict[previous_symbol]["intervals"].append(interval)
                    start_index = nucleotide[index_column_number]

            else:
                start_index = nucleotide[index_column_number]

            pseudoknot_dict[symbol]["length"] += 1

        previous_symbol = symbol

    return pseudoknot_dict


def get_pseudoknot_count(ctwuss):

    wuss_list = load_ctwuss(ctwuss)
    wuss_list_copy = deepcopy(wuss_list)
    count = 0
    previous_symbol = "0"

    for nucleotide in wuss_list_copy:
        symbol = nucleotide[symbol_column_number]

        if symbol.isalpha():
            if not previous_symbol.isalpha():
                count += 1

        previous_symbol = symbol

    return int(count/2)


# public
def save_stat(ctwuss_path, output_path, is_append=False):
    """ Saves statistics for all types of structural elements to file.
        
    Calculates and saves statistics on hairpin loops, internal loops, bulges, multifurcation loops and external loops to a specified output file.
    
    Parameters:
        ctwuss_path (str): path to the CT-modified file or loaded file
        output_path (str): path to the output file
        is_append (bool): append output to existing file or create a new file
    """

    lines = ["External_loops", "Hairpin_loops", "Internal_loops", "Bulge_loops", "Multifurcation_loops", "Stems"]
    stat = [get_external_loop_stat(ctwuss_path), get_hairpin_loop_stat(ctwuss_path),
            get_internal_loop_stat(ctwuss_path),
            get_bulge_loop_stat(ctwuss_path), get_multifurcation_loop_stat(ctwuss_path), get_stem_stat(ctwuss_path)]

    if isfile(output_path) and not is_append:
        print("This file is already exist. Delete this file and try again.")
        return -1

    output_file = open(output_path, "a")
    if is_append:
        output_file.write("\n" * 2)

    line1 = basename(ctwuss_path)
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


def dotbracket2ct(input_path, output_path):
    """ Converts dot-bracket format file to CT format file

    Parameters:
        input_path (str): path to the dot-bracket file
        output_path (str): path to the output CT file
    """

    def read_dot_bracket(path):
        file = open(path, "r")
        name = file.readline()
        letters = list(file.readline())
        symbols = list(file.readline())
        n = len(letters)
        dot_bracket_list = []

        for i in range(n):
            nucleotide = [i + 1, symbols[i], letters[i]]
            dot_bracket_list.append(nucleotide)

        return [name, dot_bracket_list[:-1]]

    name, dot_bracket_list = read_dot_bracket(input_path)
    brackets_list = []
    for nucleotide in dot_bracket_list:
        if nucleotide[1] in ["(", ")"]:
            brackets_list.append(nucleotide)

    def delete_pair(brackets_list):

        empty = [-2] * 3
        pairs = {}
        for i in range(len(brackets_list)):

            if i == 0:
                continue

            else:
                nucleotide = brackets_list[i]
                previous_nucleotide = brackets_list[i - 1]
                if nucleotide[1] == ")" and previous_nucleotide[1] == "(":
                    l = nucleotide[0]
                    r = previous_nucleotide[0]
                    pairs[l] = r
                    pairs[r] = l
                    brackets_list[i] = empty
                    brackets_list[i - 1] = empty

        while empty in brackets_list:
            brackets_list.remove(empty)

        return pairs, brackets_list

    connections_dict = {}
    while len(brackets_list) != 0:
        pairs, brackets_list = delete_pair(brackets_list)
        connections_dict.update(pairs)

    output_file = open(output_path, "w")
    length = len(dot_bracket_list)
    output_file.write(f"{length} {name[1:-1]}\n")
    n = len(str(length)) + 1
    maximum = max(list(connections_dict.keys()))
    m = len(str(maximum)) + 1
    for i in range(length):
        line = f"{i + 1:>{n}}"
        line += " " + dot_bracket_list[i][2]
        line += f"{i:>{n}}"
        line += f"{i + 2:>{n}}"
        if i + 1 in connections_dict:
            line += f"{connections_dict[i + 1]:>{m}}"
        else:
            line += f"{0:>{m}}"
        line += f"{i + 1:>{n}}"
        output_file.write(line + "\n")
    output_file.close()


# public
def get_structure_type(ctwuss, position_index):
    """ Returns the type of structural element at given 1-based position.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
            position_index (int): genome position
        
        Returns:
            element_type (string): type of the structural element
    """
    wuss_list = load_ctwuss(ctwuss)
    return wuss_list[position_index - 1][structure_column_number]


# public
def get_structure_type_full(ctwuss, position_index):
    """ Returns a type of the structural element at a given 1-based position along with other info.
        
        Parameters:
            ctwuss (str/list): path to the CT-modified file or loaded file
            position_index (int): genome position
        
        Returns:
            element_info (dict): dictionary with the following keys:
                - structure_type - type of the secondary structure element,
                - length - element length,
                - intervals - list of the start/stop positions,
                - sequences - DNA sequences at the intervals.
    """

    wuss_list = load_ctwuss(ctwuss)
    structure = wuss_list[position_index - 1][structure_column_number]

    if structure == "external_loop":
        structures_list = get_external_loop_list(wuss_list)
    elif structure == "bulge_loop":
        structures_list = get_bulge_loop_list(wuss_list)
    elif structure == "internal_loop":
        structures_list = get_internal_loop_list(wuss_list)
    elif structure == "hairpin_loop":
        structures_list = get_hairpin_loop_list(wuss_list)
    elif structure == "stem":
        structures_list = get_stem_list(wuss_list)
    elif structure == "multifurcation_loop":
        structures_list = get_multifurcation_loop_list(wuss_list)

    structure_dictionary = get_interval(structures_list, position_index)

    sequences = get_sequences(wuss_list, structure_dictionary["intervals"], position_index)

    final_dict = {"structure_type": structure, "length": structure_dictionary["length"],
                  "intervals": structure_dictionary["intervals"], "sequences": sequences}
    if structure == "hairpin_loop":
        final_dict["stem_length"] = structure_dictionary["stem_length"]

    return final_dict

