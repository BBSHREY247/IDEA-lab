class Dog:
    species = "Canine"

    def __init__(self, name, breed):
        self.name = name
        self.breed = breed

    def bark(self):
        return f"{self.name} says Woof!"

dog1 = Dog("Buddy", "Golden Retriever")
dog2 = Dog("Rex", "German Shepherd")

print(dog1.name)       
print(dog2.bark())     
print(dog1.species)    
