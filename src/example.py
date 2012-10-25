
class Rectangle:

    def __init__(self, width, height):
        self.width = width
        self.height = height
    
    
    def area(self):
        return self.width * self.height
    


class Square(Rectangle):

    def __init__(self, length):

        Rectangle.__init__(self, length, length)


    def area(self):
        return 1
