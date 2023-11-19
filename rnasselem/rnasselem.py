import sys
from math import sqrt
from operator import itemgetter
from os.path import isfile
from os.path import basename


index_column_number = 0
connection_column_number = 4
symbol_column_number = 6


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
    while(line):
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
    
    if numlines != (i-1):
        print("Error: number of actual lines " + str(i-1) + " does not equal specified in file " + str(numlines) + ".")
        sys.exit()
    
    n = numlines

    npairs = 0
    for i in range(1,n+1):
        if ct[i] != 0 and i < ct[i]:
            npairs += 1

    print("Npairs:"+str(npairs))

    rb = [-1] * 26
    ss = [":"] * n

    npairsActual = 0
    stack = []
    stackLoops = []
    stackPseudo = []
    for i in range(1,n+1):
        if ct[i] == 0:
            stack.append(i)
        elif ct[i] > i:
            stack.append(i)
        else:
            ispair = 0
            numstr = 0
            minlevel = -1
            while len(stack)>0:
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
                        ss[j-1] = '<'
                        ss[i-1] = '>'
                    elif minlevel == -2:
                        ss[j-1] = '('
                        ss[i-1] = ')'
                    elif minlevel == -3:
                        ss[j-1] = '['
                        ss[i-1] = ']'
                    elif minlevel == -4:
                        ss[j-1] = '{'
                        ss[i-1] = '}'
                    stack.append(minlevel)
                    # loops
                    while len(stackLoops) > 0:
                        j = stackLoops.pop()
                        if numstr == 0:
                            ss[j-1] = '_'
                        elif numstr == 1:
                            ss[j-1] = '-'
                        elif numstr > 1:
                            ss[j-1] = ','
                    break
                elif ct[j] == 0:
                    if ct2[j] == 0:
                        stackLoops.append(j)
                else:
                    stackPseudo.append(j)

    return(ss)


def ct2ctwuss(ctPath, outPath):

    fr = open(ctPath)
    fo = open(outPath,"w")
    # skip header
    fr.readline()
    s = fr.readline()
    i = 0
    while(s):
        s += ' '+ss[i]+'\n'
        fo.write(s)
        s = fr.readline()
        i += 1
    fr.close()
    fo.close()


def ct2wuss(ctPath, outPath):
    
    fr = open(ctPath)
    fo = open(outPath,"w")
    # skip header
    head = readline()
    headName = head.split()[len(head)-1]
    s = fr.readline()
    i = 0
    seq = ''
    while(s):
        seq += s.split()[1]
        s = fr.readline()
    
    fo.write(">"+headName+'\n')
    fo.write(seq+'\n')
    fo.write(ss+'\n')
    
    fr.close()
    fo.close()


# public
def load_ctwuss(ctwuss_path):
    """ Loads CT-modified file with WUSS annotation in additional column.
        
        Parameters:
        ctwuss_path (str):
        
        Returns:
        genome (list[list]): list of nucleotides with pairing info
    """

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


# private
def make_indexes_dict(length, start_index, end_index):

    d = {"length": length, "start_index": start_index + 1, "end_index": end_index + 1}
    return d


# public
def get_external_loop_list(ctwuss):
    """ Returns the list of external loops.
        
        Each element of the list is dictionary with the following keys:
        start_index: start position of external loop
        end_index: end position of external loop
        length: loop length
       
       Parameters:
       ctwuss (str/list): path to the CT-modified file or loaded file

       Returns:
       external_loops (list[dict]): list of dictionaries with external loop properties
    """
    
    wuss_list = load_ctwuss(ctwuss)
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
        elements = get_external_loops_list(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_hairpin_loop_list(ctwuss):
    """ Returns the list of hairpin loops.
        
        Each element of the list is dictionary with the following keys:
        start_index: start position of hairpin loop
        end_index: end position of hairpin loop
        length: loop length
        
        Parameters:
        ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
        hairpin_loops (list[dict]): list of dictionaries with hairpin loop properties
    """

    wuss_list = load_ctwuss(ctwuss)
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
        elements = get_hairpin_loops_list(ctwuss)
    
    stat = statistics(elements)
    return stat


# private
def get_bulge_and_internal_loops_list(wuss):

    wuss_list = index_minus_one(load_ctwuss(wuss))
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
def get_internal_loop_list(ctwuss):
    """ Returns the list of internal loops.
        
        Each element of the list is dictionary with the following keys:
        start_index: start position of internal loop
        end_index: end position of internal loop
        length: loop length
        
        Parameters:
        ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
        interior_loops (list[dict]): list of dictionaries with internal loop properties
    """

    interior_loops = get_bulge_and_internal_loops_list(ctwuss)[0]
    return interior_loops


# public
def get_internal_loop_stat(wuss):
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

    stat = statistics(elements)
    return stat


# public
def get_bulge_loop_list(ctwuss):
    """ Returns the list of bulges.
        
        Each element of the list is dictionary with the following keys:
        start_index: start position of bulge
        end_index: end position of bulge
        length: loop length
        
        Parameters:
        ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
        bulges (list[dict]): list of dictionaries with bulge properties
    """

    bulges = get_bulge_and_internal_loops_list(ctwuss)[1]
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

    stat = statistics(wuss)
    return stat


# public
def get_multifurcation_loop_list(ctwuss):
    """ Returns the list of multifurcation loops.
        
        Each element of the list is dictionary with the following keys:
        start_index: start position of multifurcation loop
        end_index: end position of multifurcation loop
        length: loop length
        
        Parameters:
        ctwuss (str/list): path to the CT-modified file or loaded file
        
        Returns:
        multifurcation_loops (list[dict]): list of dictionaries with multifurcation loop properties
    """
    
    wuss_list = load_ctwuss(ctwuss)
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
        elements = get_multifurcation_loops_list(ctwuss)

    stat = statistics(elements)
    return stat


# public
def get_stem_list(ctwuss):
    """ Returns the list of stems.
        
        Each element of the list is dictionary with the following keys:
        start_index: start position of stem
        end_index: end position of stem
        length: loop length
        
        Parameters:
        ctwuss (str): path to the CT-modified file or loaded file
        
        Returns:
        stem_list (list[dict]): list of dictionaries with stems properties
    """

    wuss_list = load_ctwuss(ctwuss)
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

    stat = statistics(elements)
    stat["total_length"] *= 2

    return stat


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
    stat = [get_external_loop_stat(ctwuss_path), get_hairpin_loop_stat(ctwuss_path), get_internal_loop_stat(ctwuss_path),
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
