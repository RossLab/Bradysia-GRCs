from sys import stderr
from sys import stdout
from sys import argv

depth_filename = argv[1]
# depth_filename = 'data/L-X-A-kmers/mapping/illumina_inspected_scfs.depth'
window_size = int(argv[2])
# window_size = 5000

# scf from to L X A

with open(depth_filename) as depth_file:
    previous_scf = ""
    window_Ls = 0
    window_Xs = 0
    window_As = 0
    window = 0
    last_pos = 0
    for line in depth_file:
        scf, pos, L, X, A = line.rstrip('\n').split('\t')
        if previous_scf != scf and previous_scf != "":
            w_from = (window * window_size) + 1
            w_to = int(last_pos)
            stdout.write(("{}\t"+ "{:.0f}\t" * 4 + "{:.0f}\n").format(scf, w_from, w_to, window_Ls, window_Xs, window_As))
            previous_scf = scf
            window = 0
        if ((window + 1) * window_size) < int(pos):
            # write_down_result
            w_from = (window * window_size) + 1
            w_to = ((window + 1) * window_size)
            stdout.write(("{}\t"+ "{:.0f}\t" * 4 + "{:.0f}\n").format(scf, w_from, w_to, (window_Ls / 27), (window_Xs / 27), (window_As / 27)))
            window_Ls = 0
            window_Xs = 0
            window_As = 0
            window += 1
        else:
            window_Ls += int(L)
            window_Xs += int(X)
            window_As += int(A)
            last_pos = pos
            previous_scf = scf

# dynamic block output
# with open(depth_filename) as depth_file:
#     previous_scf = ""
#     prev_pos = 1
#     prev_state = 'NA'
#     for line in depth_file:
#         scf, pos, L, X, A = line.rstrip('\n').split('\t')
#         previous_scf != scf:
#             # write_down_result
#             previous_scf = scf
#             prev_pos = int(pos) + 1
#             continue
#         state = (int(L) > 0) * 'L' + (int(X) > 0) * 'X' + (int(A) > 0) * 'A'
#         prev_state != state and prev_state != '':
#             # write_down_result
#             prev_state = 'NA'
