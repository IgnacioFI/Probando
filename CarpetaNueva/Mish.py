class Persona:

    def __init__(self, nombre: str):
        self.nombre = nombre
    def __str__(self):
        return self.nombre


persona = Persona("Ignacio")
nombres = ["Valentin", "Barbara", "Patricia"]
grupo = []
grupo.append(persona)
print(grupo, grupo[0])
for name in nombres:
    persona = Persona(name)
    grupo.append(persona)
    print(persona)
print(grupo)
    
