import os, sys, glob
from scipy.spatial import distance

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

def eucl_dist(x, y):
    dst = distance.euclidean(x, y)
    return dst

def get_simulation_dist(instance_dir, replication, emews_root):
    """
    @return distance value between fname and "data/original_physiboss_timeseries", or -2 if file doesn't exist, or
    -1 if run terminated prematurely.
    """

    output = '-2'
    fname = '{}/output/metrics.txt'.format(instance_dir)
    if os.path.exists(fname):
        file_lines = []
        with open(fname) as f_in:
            output = '-1'
            tmp = f_in.readlines()
            file_lines = [x.strip() for x in tmp]
        # file_lines.reverse() 
        check = file_lines[-1].split("\t")
        if len(check) > 1:
            tumor_cells = []
            death_cells = []
            necrosis_cells = []
            for i in range(len(file_lines)):
                items = file_lines[i].split("\t")
                tumor_cells.append(int(items[1]))
                death_cells.append(int(items[2]))
                necrosis_cells.append(int(items[3]))

            # Find and parse corresponding original csv
            for i, f in enumerate(sorted(glob.glob(emews_root+'/data/original_physiboss_timeseries/*.csv'))):
                if i == int(replication):
                    csv_file = f
            fh = open(csv_file)
            tmp = fh.readlines()
            lines = [x.strip() for x in tmp]
            alive = []
            apoptotic = []
            necrotic = []
            for i in range(len(lines)):
                if i>0:
                    data = lines[i].split('\t')
                    alive.append(float(data[2]))
                    apoptotic.append(float(data[3]))
                    necrotic.append(float(data[4]))
            lmin = min(len(alive),len(tumor_cells))
            
            if lmin > 0:
                
                alive = alive[0:lmin]
                apoptotic = apoptotic[0:lmin]
                necrotic = necrotic[0:lmin]

                tumor_cells = tumor_cells[0:lmin]
                death_cells = death_cells[0:lmin]
                necrosis_cells = necrosis_cells[0:lmin]

                output = eucl_dist(alive, tumor_cells)
                output += eucl_dist(apoptotic, death_cells)
                output += eucl_dist(necrotic, necrosis_cells)

    return output
