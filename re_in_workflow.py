# RE_IN Workflow
# We calculate statistic-based inferences on multiple matrices 
# to return a single list of interactions which RE_IN can interpret
# to build a gene regulatory boolean network. 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

conditions = ['2iLif', 'CiTC', 'LB']

### CALCULATE INFERENCE
# input : string matrix_name - ex:'CiTC' or '2iLif'
# input : string inference_type -
#   'pearson_corr' for Pearson correlation coefficients
#   'bayes' for Bayesian-based inference
#   'reg' for regression-based inference 
# input : threshold - retain the metric above that value (default 0.3)
# output : list of interactions for chosen genes in genelist

def calculate_inference(matrix_name, inference_type, threshold=0.3):

    ### DATA IMPORT
    filename = 'scRecover+scImpute_' + matrix_name + '.csv'
    df = pd.read_csv(filename, index_col=0)
    
    # genelist.csv contains gene_name + update_function + expression in each condition
    genedata = pd.read_csv('genelist.csv', header=0)
    genelist = genedata['gene_name'].values.tolist()
    
    # subset df from genelist
    df = df.loc[genelist]
    df = df.transpose()
    
    
    ### STATISTIC-BASED INFERENCE
    # user chooses statistic-based inference type
    
    ### PEARSON CORRELATION MATRIX
    if(inference_type == 'pearson_corr'):
        
        #correlation of matrix
        df_corr = df.corr(method='pearson')   
        filename = 'results/pearson_corr_' + matrix_name + '.csv'
        df_corr.to_csv(filename)
        
        #generate heatmap 
        f = plt.figure(figsize=(25, 20))
        plt.matshow(df_corr, fignum=f.number, cmap='PiYG')
        plt.xticks(range(df_corr.shape[1]), df_corr.columns, fontsize=9,rotation=90)
        plt.yticks(range(df_corr.shape[1]), df_corr.columns, fontsize=9)
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=17)
        
        plt.title('Pearson correlation', fontsize=20);
        
        #graph representation
        filename = 'results/pearson_heatmap_' + matrix_name + '.png'
        plt.savefig(filename,bbox_inches='tight')
    
    ### BAYESIAN-BASED INFERENCE
    if(inference_type == 'bayes'):
        print('bayesian-based inference')
    ### REGRESSION-BASED INFERENCE
    if(inference_type == 'reg'):
        print('regression-based inference')
    
    
    
    ### EXTRACT LIST OF INTERACTIONS
    
    pos = np.where(df_corr >= threshold)
    pos = [(df_corr.index[x], df_corr.columns[y], 'positive', 'optional') for x, y in zip(*pos)
                                            if x != y and x < y]
    filename = 'results/positive_corr_' + matrix_name + '.csv'
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(pos)
        
    interac_list = pd.DataFrame(data=pos, columns=['gene_source', 'gene_target', 'interaction_type', 'interaction_definite'])    
    
    return interac_list

interac_2iLif = calculate_inference(matrix_name='2iLif', inference_type='pearson_corr', threshold=0.3)
interac_CiTC = calculate_inference(matrix_name='CiTC', inference_type='pearson_corr', threshold=0.3)
interac_LB = calculate_inference(matrix_name='LB', inference_type='pearson_corr', threshold=0.3)



### MERGE ALL LISTS OF INTERACTIONS

interac_list_merge = [interac_2iLif, interac_CiTC, interac_LB]
interac_list_merge = pd.concat(interac_list_merge)
interac_list_merge.index = range(len(interac_list_merge))
interac_list_merge = interac_list_merge.drop_duplicates(['gene_source','gene_target'],keep= 'first')



### UPDATE LIST OF INTERACTIONS

definite = pd.read_csv('definite_interactions.csv')


negative = pd.read_csv('negative_interactions.csv')

remove = pd.read_csv('remove_interactions.csv')



def change_interaction_type(interac_list):
    
    try:
        lists=interac_list.values.tolist()
        newl=negative.values.tolist()
        
        for i in lists:
            tmp=[i[0],i[1]]
            if tmp in newl:
                i[2]='negative'
        

        interac_lists=lists

    except pd.errors.EmptyDataError:
        return 0


def remove_optional(interac_list):

    try:
        lists=interac_list.values.tolist()
        newl=definite.values.tolist()
    
        for i in lists:
            tmp=[i[0],i[1]]
            if tmp in newl:
                i.remove(i[3])
        
     
        interac_list=lists
    except pd.errors.EmptyDataError:
        return 0
            

def remove_interaction(interac_list):
    lists=interac_list.values.tolist()
    newl=remove.values.tolist()
    
    try:
        for i in lists:
            tmp=[i[0],i[1]]
            if tmp in newl:
                lists.remove(i)

    except pd.errors.EmptyDataError:
        return 0
    
            

    
def write_positive_csv(interac_list):
    
        lists=interac_list.values.tolist()
    
        for line in lists:
                line[len(line)-1]=str(line[len(line)-1])+';'
    
        with open("results/positive_corr.csv", "w", newline="") as g:
            writer = csv.writer(g)
            writer.writerows(lists) 
    
    
change_interaction_type(interac_list_merge)
remove_optional(interac_list_merge)
remove_interaction(interac_list_merge)
write_positive_csv(interac_list_merge)

### RE:IN INTERACTION FILE

pos = pd.read_csv('results/positive_corr.csv')
pos= pos.values.tolist()

gene_name=[]
reg_cond=[]
cond1=[]
cond2=[]
cond3=[]
txt_list=[]

with open('genelist.csv','r') as csv_file:
    contents = [x.strip() for x in csv_file.readlines()[1:]]
  

for line in contents:
    data = line.split(',')
    print(data[0])
    gene_name.append(data[0])
    
    reg_cond.append(data[1])
    cond1.append(data[2])
    cond2.append(data[3])
    cond3.append(data[4])
    

for row in pos:
    print(row)
    if (row[0] in gene_name)& (row[1] in gene_name):

        txt_list.append(" ".join(row))
print(txt_list)
    
txt_file="results/interactions.txt"

with open(txt_file, "w") as my_output_file:
        my_output_file.write("// Settings\ndirective updates sync;\ndirective length 20;\ndirective regulation noThresholds;\n\n// Components\n")
        
        for i in range (len(gene_name)):
            my_output_file.write(gene_name[i]+"[]"+reg_cond[i]+";\n")
        
        my_output_file.write("\n")
        for line in txt_list:
            my_output_file.write(line)
            my_output_file.write("\n")



### RE:IN OBSERVATION FILE

text_file="results/observations.txt"
with open(text_file, "w") as my_output_file:
    my_output_file.write("//Observation predicates\n\n$Expression2iLif\n:=\n{\n")
    
    for i in range (len(gene_name)-1):
        my_output_file.write(gene_name[i]+" = "+cond1[i]+" and\n")
    my_output_file.write(gene_name[len(gene_name)-1]+" = "+cond1[len(gene_name)-1]+"\n")
    
    my_output_file.write("};\n\n\n$ExpressionCiTC\n:=\n{\n")
    
    for i in range (len(gene_name)-1):
        my_output_file.write(gene_name[i]+" = "+cond2[i]+" and\n")
    my_output_file.write(gene_name[len(gene_name)-1]+" = "+cond2[len(gene_name)-1]+"\n")
   
    my_output_file.write("};\n\n\n$ExpressionLB\n:=\n{\n")
    for i in range (len(gene_name)-1):
        my_output_file.write(gene_name[i]+" = "+cond3[i]+" and\n")
    my_output_file.write(gene_name[len(gene_name)-1]+" = "+cond3[len(gene_name)-1]+"\n")
    
    my_output_file.write("};\n\n\n// Observations\n#Experiment2iLif[0] |= $Expression2iLif;\n#ExperimentCiTC[0] |= $ExpressionCiTC;\n#ExperimentLB[0] |= $ExpressionLB;")
                             
                                                                             







