import extinction

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']
    Filter = ['V']

    # Package the inputs for extinction.py
    inputs = (Dataset[0],sets,Filter[0])
    
    # Run the reduce script
    numstars, fileout = extinction.extinction(*inputs)

if __name__ == "__main__":
    main()