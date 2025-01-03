import register

def main():

    # Define the 
    Dataset = ["ROMO241004"]
    sets = ['1st']
    Filter = ['V']

    # Package the inputs for register.py
    inputs = (Dataset[0],sets,Filter[0]) 
    
    # Run the reduce script
    cropped_fn, failed_fn = register.matchstars(*inputs)

if __name__ == "__main__":
    main()