import coordinates

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']

    # Package the inputs for pointing.py
    inputs = (Dataset[0],sets)
    
    # Run the coordinates script
    coordinates.galactic_ecliptic_coords(*inputs)

if __name__ == "__main__":
    main()