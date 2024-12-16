import pointing

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']

    # Package the inputs for pointing.py
    inputs = (Dataset[0],sets)
    
    # Run the pointing script
    pointing.pointing_err(*inputs)

if __name__ == "__main__":
    main()