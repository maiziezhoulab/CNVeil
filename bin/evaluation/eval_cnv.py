import os 
import pandas as pd
import numpy as np
from collections import defaultdict
def find_boundaries(numbers):
    boundaries = []
    for i in range(1, len(numbers)):
        if numbers[i] != numbers[i-1]:
            boundaries.append(i)
    return boundaries

def extract_chrom_boundary(df):
    chr_list =[ x.split(':')[0] for x in df.columns.tolist()]
    chr_bound_list  = find_boundaries(chr_list)
    return chr_bound_list

def extract_cnv_boundary(cnv_list,chr_bound_list):
    raw_bound_list =  find_boundaries(cnv_list)
    # print(raw_bound_list)
    final_bound_list = sorted(list(set(raw_bound_list) - set(chr_bound_list)))
    # print(final_bound_list)
    return final_bound_list 

def format_bnd(final_bound_list, bin_list, cnv_list, cell_name):

    bnd_list =  [
        (
        bin_list[idx-1].split(':')[0], 
        int(bin_list[idx-1].split(':')[1].split('-')[1]), 
        int(cnv_list[idx-1]), 
        int(cnv_list[idx]),
        cell_name
        )
        for idx in final_bound_list
        ]

    return bnd_list

def write_bnds(bnd_list,outfile,file_type ):
    try:
        assert file_type in ['comp','bench','bed']
    except:
        print("file type can only be comp, bench or bed")
        exit()
    with open(outfile,'w') as f:
        if file_type == 'comp':
            line = f"chrom\tpos\tleft_cn\tright_cn\tcell_name\n"
            f.write(line)
        elif file_type == 'bench':
            line = f"chrom\tpos\tleft_cn\tright_cn\tnode\n"
            f.write(line)
        for bnd in bnd_list:
            chrom,pos,left_cn, right_cn, cell_name = bnd
            line = f"{chrom}\t{pos}\t{left_cn}\t{right_cn}\t{cell_name}\n"
            f.write(line)



def load_comp(input_path, output_dir):
    os.system("mkdir -p " + output_dir)
    prefix  = os.path.basename(input_path).split('.')[0]
    outfile = output_dir + f'/comp_bnd.csv'
    df = pd.read_csv(input_path, index_col = 0)
    chr_bound_list =  extract_chrom_boundary(df)
    bin_list =  df.columns.tolist()
    bnd_all_cells = []
    for i in range(df.shape[0]):
        cell_name = df.index[i]
        # print(cell_name)
        if 'normal' not in cell_name:
            cnv_list = list(df.iloc[i,:].values)
            final_bnd_list = extract_cnv_boundary(cnv_list,chr_bound_list)
            # print(cnv_list,len(final_bnd_list))
            

            bnd_all_cells.extend(format_bnd(final_bnd_list, bin_list, cnv_list, cell_name))

    write_bnds(bnd_all_cells, outfile,'comp')
    # print('xxx',df.shape,len(bnd_all_cells))
    return bnd_all_cells, outfile

def load_fai(fai_path):

    cnt = 0
    dc = {}
    with open(fai_path,'r') as f:
        for line in f:
            if cnt>= 22:
                break
            data = line.split()
            chrom , length = data[0], data[1]
            dc[chrom] = length

            cnt+=1
    return dc


def fill_in_gap(bedfile, ref_type,outdir, int_chrom_set, leafnodes):


    # set fai path 
    try:
        assert ref_type in ['19','38']
    except:
        print("ref type can only be 19 or 38")
        exit()
    if ref_type =='19':
        fai_path = "/data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa.fai"
    else:
        fai_path = "/data/maiziezhou_lab/Softwares/refdata-GRCh38-2.1.0/fasta/genome.fa.fai"

    dc_chrom = load_fai(fai_path)
    leafnodes = sorted(list(leafnodes))
    dc_list = [defaultdict(list) for i in range(len(leafnodes))]
    with open(bedfile,'r') as f:
        for line in f:
            fields = line[:-1].split()
            if (fields[0] in int_chrom_set) & ( fields[-1] in leafnodes):
                leaf_idx = leafnodes.index(fields[-1])
                dc_list[leaf_idx][fields[0]].append(fields)
    final_list = []
    gap_cnt = 0
    for leaf_idx in range(len(dc_list)):
        node = leafnodes[leaf_idx]
        dc = dc_list[leaf_idx]
        for chrom in int_chrom_set:
            gapped_list = dc[chrom]
            # print(node, chrom, gapped_list)
            if gapped_list:
                gapless_list = []
                if gapped_list[0][1]!='0':
                    gapless_list.append((chrom, 0,gapped_list[0][1],2, node))
                gapless_list.append(gapped_list[0])

                for i in range(1, len(gapped_list)):
                    prev_end = gapped_list[i-1][2]
                    cur_start = gapped_list[i][1]
                    if cur_start != prev_end:

                        gapless_list.append((chrom, prev_end, cur_start, 2, node))
                    gapless_list.append(gapped_list[i])
                
                chrom_len = dc_chrom[chrom]
                if gapped_list[-1][2]!=str(chrom_len):
                    gapless_list.append((chrom, gapped_list[-1][2], chrom_len, 2, node))

                gap_cnt+= len(gapless_list)-len(gapped_list)

                final_list.extend(gapless_list)
    outfile = outdir + '/bench_no_gap.bed'
    write_bnds(final_list ,outfile,'bed' )
    print(f"filled {gap_cnt} gaps in bench bed")
    return final_list


def load_bench(bench_dir, outdir,ref_type):
    bedfile = bench_dir+'/gt_all.bed'
    leafnodefile = bench_dir+'/all_leaf_node'
    os.system("mkdir -p "+ outdir)
    with open(leafnodefile,'r') as f:
        leafnodes = list(f.read().split('\n')[:-1])
    int_chrom_set = set(['chr'+str(i) for i in range(1,23)])


    ### fill in non-diploid gaps in bed file

    bed_list = fill_in_gap(bedfile, ref_type ,outdir, int_chrom_set, leafnodes)



    bnd_list = []
    for i in range(1,len(bed_list)):
        prev_node = bed_list[i-1][-1]
        cur_node = bed_list[i][-1]
        prev_chrom = bed_list[i-1][0]
        cur_chrom = bed_list[i][0]
        prev_cn = int(bed_list[i-1][3])
        cur_cn = int(bed_list[i][3])
        if (prev_node == cur_node) & (prev_chrom == cur_chrom) & ( prev_cn != cur_cn):
            prev_end = bed_list[i-1][2]
            cur_start = bed_list[i][1]
            try:
                assert prev_end == cur_start
            except:
                print(f"node{cur_node}, gap {cur_chrom}:{prev_end}-{cur_start} is not filled, should be filled with copy number 2")
                exit()
            
            cur_pos = int(bed_list[i][1])
            bnd_list.append((cur_chrom, int(cur_pos), int(prev_cn), int(cur_cn), 'node'+cur_node))

    write_bnds(bnd_list, outdir+'/bench_bnd.csv', 'bench')

    return bnd_list,outdir+'/bench_bnd.csv'


def reconstruct_comp_list(comp_list):

    clusters = [[comp_list[0]]]

    for i in range(1, len(comp_list)):
        cur_cname = comp_list[i][-1]
        prev_cname = comp_list[i-1][-1]
        if cur_cname == prev_cname:
            clusters[-1].append(comp_list[i])
        else:
            clusters.append([comp_list[i]])

    return clusters


def reconstruct_bench_list(bench_list):
    dc = defaultdict(list)

    for bed in bench_list:
        node = bed[-1]
        dc[node].append(bed)

    return dc 

def compare_one_bench_all_comp(bench_bnd, comp_bnd_list,used_comp_id_set, max_dist_thresh, match_cn):
    chrom_bench, pos_bench, lcn_bench, rcn_bench , node_bench = bench_bnd 
    cand_comp_ids  = []
    dist_list = []

    for i in range(len(comp_bnd_list)):
        if i not in used_comp_id_set:
            comp_bnd = comp_bnd_list[i]
            chrom_comp, pos_comp, lcn_comp, rcn_comp, cell_name = comp_bnd
        
            if chrom_bench!= chrom_comp:
                continue

            dist = abs(pos_bench-pos_comp)
            if dist > max_dist_thresh:
                continue
            if match_cn:
                if (lcn_comp, rcn_comp)  != (lcn_bench, rcn_bench):
                    continue
            cand_comp_ids.append(i)
            dist_list.append(i)
    if cand_comp_ids:
        min_d  = min(dist_list)
        min_id = dist_list.index(min_d)
        return cand_comp_ids[min_id]
    else:
        return -1

        




        

def compare_two_bnd_list(comp_bnd_list, bench_bnd_list, max_dist_thresh, match_cn):
    used_comp_ids = []
    tp = 0
    fn = 0
    matching_info = [ ]
    for i in range(len(bench_bnd_list)):
        bench_bnd = bench_bnd_list[i]
        match_comp_id = compare_one_bench_all_comp(bench_bnd, comp_bnd_list, used_comp_ids, max_dist_thresh, match_cn)
        if match_comp_id>0:
            tp+=1
            used_comp_ids.append(match_comp_id)
            matching_info.append( bench_bnd + comp_bnd_list[match_comp_id])
        else:
            fn+=1

    assert (tp+fn) == len(bench_bnd_list)

    fp = len(comp_bnd_list )  - tp 
    return tp,fp,fn,matching_info

def write_matching_info(matching_list, outdir, match_cn, max_dist_thresh):
    if match_cn:
        outfile = outdir+f'/matching_info_{max_dist_thresh}_matchcn.txt'
    else:
        outfile = outdir+f'/matching_info_{max_dist_thresh}.txt'

    cols = ['chrom_bench','pos_bench','lcn_bench','rcn_bench','node_bench',
                    'chrom_comp','pos_comp','lcn_comp','rcn_comp','cell_name_comp']
    if matching_list:

        df = pd.DataFrame(matching_list)
        # print(df)
        df.columns = cols
        df.to_csv(outfile, index = False)
    else:
        with open(outfile,'w') as f:
            f.write(','.join(cols)+'\n')


    return 


def write_stats(stats, outdir,match_cn, max_dist_thresh):
    if match_cn:
        outfile = outdir+f'/stats_by_cell_{max_dist_thresh}_matchcn.csv'
    else:
        outfile = outdir+f'/stats_by_cell_{max_dist_thresh}.csv'
    df = pd.DataFrame(stats)
    df.columns = ['tp','fp','fn','node','cell']
    df['recall'] = np.round(df['tp']/ (df['tp']+df['fn']),4)
    df['precision'] = np.round(df['tp']/ (df['tp']+df['fp']),4)
    avg_recall = np.round(df['recall'].mean(),4)
    avg_precision = np.round(df['precision'].mean(),4)
    log = f'match copy number = {match_cn} max_dist_thresh = {max_dist_thresh}: avg recall {avg_recall}, avg precision {avg_precision}\n'
    print(log)
    df.to_csv(outfile, index = False)
    logfile = outfile.replace('.csv','.log')
    with open(logfile,'w') as f:
        f.write(log)
    return 
        
def eval_cnv(bench_dir, comp_file, outdir , max_dist_thresh, ref_type, match_cn):

    comp_list , comp_outfile = load_comp(comp_file, outdir)

    bench_list , bench_outfile = load_bench(bench_dir, outdir, ref_type)
    # print(comp_list[0])
    # print(bench_list[0])


    comp_list_by_cell = reconstruct_comp_list(comp_list)
    dc_bench = reconstruct_bench_list(bench_list)
    matching_info_all_cell = []
    sta_all_cell = []
    for comp_list_one_cell in comp_list_by_cell:
        node, cell = comp_list_one_cell[0][-1].split('_')
        bench_list_one_node = dc_bench[node]
        tp,fp,fn,matching_info = compare_two_bnd_list(comp_bnd_list= comp_list_one_cell,
                                                      bench_bnd_list= bench_list_one_node,
                                                      max_dist_thresh= max_dist_thresh,
                                                      match_cn= match_cn)
        
        sta_all_cell.append([tp,fp,fn, node, cell])
        matching_info_all_cell.extend(matching_info)


    write_matching_info(matching_info_all_cell, outdir, match_cn, max_dist_thresh)

    write_stats(sta_all_cell, outdir, match_cn, max_dist_thresh)

    return 





        
        




if __name__ == '__main__':

    import argparse
    from argparse import ArgumentParser

    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--comp_file','-c')
    parser.add_argument('--bench_dir','-b')
    parser.add_argument('--output_dir','-o')
    parser.add_argument('--max_dist_thresh',  '-d', default = 2000000)
    parser.add_argument('--ref_type',  '-r', default = '19', choices= ['19','38'])
    # parser.add_argument('--match_cn','-m', action='store_true')
    args = parser.parse_args()
    comp_file = args.comp_file
    bench_dir = args.bench_dir
    output_dir = args.output_dir
    max_dist_thresh = args.max_dist_thresh
    ref_type = args.ref_type
    # match_cn = args.match_cn

    if not os.path.exists(comp_file):
        print(f"{comp_file} not exists.")
        exit()


    eval_cnv(bench_dir, comp_file, output_dir , max_dist_thresh, ref_type, match_cn = True)
    eval_cnv(bench_dir, comp_file, output_dir , max_dist_thresh, ref_type, match_cn = False)





