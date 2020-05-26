import os, sys, glob

def get_tumor_cell_count(instance_dir):
    """
    @return tumor cell count value from fname or -2, if file doesn't exist, or
    -1 if run terminated prematurely.
    """
    tumor_cell_count = '-2'
    fname = '{}/output/metrics.txt'.format(instance_dir)
    if os.path.exists(fname):
        file_lines = []
        with open(fname) as f_in:
            tumor_cell_count = '-1'
            file_lines.append(f_in.readlines()[-1].strip())
        file_lines.reverse() 
        items = file_lines[0].split("\t")
        if len(items) > 1:
            tumor_cell_count = items[1]

    return tumor_cell_count

def get_custom_cell_count(instance_dir):
    """
    @return tumor cell count value from fname or -2, if file doesn't exist, or
    -1 if run terminated prematurely.
    """
    output = '-2'
    fname = '{}/output/metrics.txt'.format(instance_dir)
    if os.path.exists(fname):
        file_lines = []
        with open(fname) as f_in:
            output = '-1'
            file_lines.append(f_in.readlines()[-1].strip())
        file_lines.reverse() 
        items = file_lines[0].split("\t")
        if len(items) > 1:
            tumor_cell_count = int(items[1])
            death_cell_count = int(items[2])
            necrosis_cell_count = int(items[3])
            total_cells = tumor_cell_count + death_cell_count + necrosis_cell_count

            tumor_percent = tumor_cell_count * 100 / total_cells
            necrosis_percent = necrosis_cell_count * 100 / total_cells
            output = tumor_percent - necrosis_percent

    return output
