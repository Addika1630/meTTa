from hyperon import MeTTa, SymbolAtom, ExpressionAtom, GroundedAtom
import os
import glob

metta = MeTTa()
metta.run(f"!(bind! &space (new-space))")

def load_dataset(path: str) -> None:
    if not os.path.exists(path):
        raise ValueError(f"Dataset path '{path}' does not exist.")
    paths = glob.glob(os.path.join(path, "**/*.metta"), recursive=True)
    if not paths:
        raise ValueError(f"No .metta files found in dataset path '{path}'.")
    for path in paths:
        print(f"Start loading dataset from '{path}'...")
        try:
            metta.run(f'''
                !(load-ascii &space {path})
                ''')
        except Exception as e:
            print(f"Error loading dataset from '{path}': {e}")
    print(f"Finished loading {len(paths)} datasets.")

# Example usage:
try:
    dataset = load_dataset("./Data")
   
except Exception as e:
    print(f"An error occurred: {e}")

# 2 Points
def get_transcript(node):
    
    """
    Fetches transcripts associated with a specified gene identifier.

    Parameters:
    node (list): A list containing the gene identifier as its first element.

    The function constructs a MeTTa query to find all transcripts 
    that are transcribed from the given gene. The query matches 
    expressions of the form (transcribed_to gene transcript). 

    Returns:
    gene: The result of the MeTTa query, containing the associated transcripts.
    """

    gene = metta.run(f''' 
        !(match &space 
        (,(transcribed_to ({node[0]}) $transcript)) 
           (,(transcribed_to ({node[0]}) $transcript))) 
         
        ''') 
    return gene       

#2 Points
def get_protein(node):
    """
    Fetches proteins associated with a specified transcript identifier.

    Parameters:
    node (list): A list containing the gene identifier as its first element.

    The function constructs a MeTTa query to find all proteins 
    that are translated from the given transcript, following the 
    relationship defined by (translates_to transcript protein).

    Returns:
    protein: The result of the MeTTa query, containing the associated proteins.
    """

    protein = metta.run(f''' 
        !(match &space 
        (, (transcribed_to ({node[0]}) $transcript) 
        (translates_to $transcript $p) 
            ) 
           (, (translates_to $transcript $p))) 
         
        ''') 
    return protein 

#6 Points
def metta_seralizer(metta_result):
    """
    Converts the output of a MeTTa query into a structured format for easy serialization.

    Parameters:
    metta_result (list): The result of a MeTTa query, expected to contain ExpressionAtoms.

    The function iterates through the provided MeTTa result, extracting edges, sources, and 
    targets from ExpressionAtoms. It formats these values as strings, removing any surrounding 
    parentheses, and appends them to a list of dictionaries.

    Returns:
    list: A list of dictionaries, each containing 'edge', 'source', and 'target' keys 
    representing the relationships extracted from the MeTTa result.
    """

    metta_result = metta_result[0]  
    result = []  

    for expr in metta_result:
        if isinstance(expr, ExpressionAtom):
            for child in expr.get_children():
                if isinstance(child, ExpressionAtom):
                    edge = child.get_children()[0]  
                    source = child.get_children()[1]  
                    target = child.get_children()[2]  
                    
                    edge_value = str(edge) 
                    source_value = str(source)  
                    target_value = str(target)  

                    result.append({
                        'edge': edge_value.strip("()"),  
                        'source': source_value.strip("()"),  
                        'target': target_value.strip("()")  
                    })

    return result



#1
transcript_result= (get_transcript(['gene ENSG00000166913']))
print(transcript_result) 
"""
Expected Output Format::
# [[(, (transcribed_to (gene ENSG00000166913) (transcript ENST00000372839))), (, (transcribed_to (gene ENSG00000166913) (transcript ENST00000353703)))]]
""" 

#2
protein_result= (get_protein(['gene ENSG00000166913']))
print(protein_result) 
"""
Expected Output Format::
# [[(, (translates_to (transcript ENST00000353703) (protein P31946))), (, (translates_to (transcript ENST00000372839) (protein P31946)))]]
"""

#3
parsed_result = metta_seralizer(transcript_result)
print(parsed_result) 
"""
Expected Output Format:
[
    {'edge': 'transcribed_to', 'source': 'gene ENSG00000175793', 'target': 'transcript ENST00000339276'}
]
"""

