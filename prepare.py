import argparse
import os


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-source', type=str, required=True)
    parser.add_argument('-score', type=float,default=0.8)
    parser.add_argument('-t', type=int,default=2)
    parser.add_argument('-kl', type=int,default=5)
    parser.add_argument('-ku', type=int,default=8)
    parser.add_argument('-more', type=int, default=0)
    parser.add_argument('-template',type=str,default='homo')
    parser.add_argument('-predfull_path',type=str,required=True)
    parser.add_argument('-rapsearch_path',type=str,required=True)
    parser.add_argument('-prerapsearch_path', type=str, required=True)
    parser.add_argument('-msslash_path',type=str,required=True)
    args = parser.parse_args()
    return args




if __name__ == '__main__':
    args = get_args()
    predfull_path = args.predfull_path
    rapsearch_path = args.rapsearch_path
    prerapsearch_path = args.prerapsearch_path
    msslash_path = args.msslash_path
    froot = f'{args.source}_{args.kl}-{args.ku}mer_{args.score}_{args.t}'
    os.system(f'python I_generateInputReads.py -source {args.source} -score {args.score} -t {args.t} -kl {args.kl} -ku {args.ku} -more {args.more} -predfull_path {predfull_path} -msSlash_path {msslash_path}')
    os.system(f'python II_assembleFromReads.py -froot {froot}')
    os.system(f'python III_sortOutputs.py -froot {froot}')
    os.system(f'python IV_matchRegion.py -froot {froot} -template {args.template} -rapsearch_path {rapsearch_path} -prerapsearch_path {prerapsearch_path}')
    os.system(f'python V_FindTemplateGroups.py -froot {froot} -source {args.source} -rapsearch_path {rapsearch_path} -prerapsearch_path {prerapsearch_path}')
    os.system(f'python concatenate.py -froot {froot} -chain light -rapsearch_path {rapsearch_path} -prerapsearch_path {prerapsearch_path}')
    os.system(f'python concatenate.py -froot {froot} -chain heavy -rapsearch_path {rapsearch_path} -prerapsearch_path {prerapsearch_path}')