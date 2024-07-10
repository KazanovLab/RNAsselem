import pytest
from rnasselem import save_stat
import os

output_path = "test_result1.txt"
path_1 = r"tests/SARS-CoV-2_34168138_tests.ctwuss"
path_2 = r"tests/HCV_H77_1725_test.ctwuss"
path_3 = r"tests/DENV_test.ctwuss"
result_1 = [6, 18, 24, 22, 3, 21]
result_2 = [7, 20, 28, 25, 6, 26]
result_3 = [13, 31, 18, 12, 8, 39]


@pytest.mark.parametrize("input_file, results_list", [
    (path_1, result_1), (path_2, result_2), (path_3, result_3)
])
def test_stat(input_file, results_list):
    count = []
    save_stat(input_file, output_path)
    output_file = open(output_path, "r")
    output_file.readline()
    output_file.readline()

    for i in range(6):
        count.append(output_file.readline().split("\t")[1])
        count = [int(i) for i in count]

    output_file.close()
    os.remove(output_path)
    assert count == results_list
