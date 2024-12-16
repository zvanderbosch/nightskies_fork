import reduce

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']
    Filterset = {
        'V': 'ML3Flat_V_2023_1.fit',
        'B': 'ML3Flat_B_2023_1.fit'
    }
    Curve = ['ML3new']

    # Package the inputs
    inputs = (Dataset[0],sets,Filterset['V'],Curve[0])
    
    # Run the reduce script
    reduce.reducev(*inputs)

if __name__ == "__main__":
    main()