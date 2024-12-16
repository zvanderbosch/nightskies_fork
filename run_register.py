import register

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']
    Filterset = {
        'V': 'ML3Flat_V_2023_1.fit',
        'B': 'ML3Flat_B_2023_1.fit'
    }
    Filter = ['V']
    Curve = ['ML3new']

    # Package the inputs
    inputs = (Dataset[0],sets,Filter[0])                     # for register.py
    
    # Run the reduce script
    cropped_fn, failed_fn = register.matchstars(*inputs)

if __name__ == "__main__":
    main()