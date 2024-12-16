import pointing

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']

    # Package the inputs
    inputs = (Dataset[0],sets) # for pointing.py
    
    # Run the pointing script
    pointing.matchstars(*inputs)

if __name__ == "__main__":
    main()