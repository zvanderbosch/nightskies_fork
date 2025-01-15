import fullmosaic

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']
    Filter = ['V']

    # Package the inputs for extinction.py
    inputs = (Dataset[0],sets,Filter[0])
    
    # Run the coordinates script
    fullmosaic.mosaic(*inputs)

if __name__ == "__main__":
    main()