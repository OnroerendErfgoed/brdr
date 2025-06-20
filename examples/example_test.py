from brdr.enums import OpenDomainStrategy

value = OpenDomainStrategy["AS_IS"
        ].value

name = OpenDomainStrategy["AS_IS"
        ].name

print (value)
print (name)
print (OpenDomainStrategy(0).name)
print (OpenDomainStrategy.AS_IS.name)
print (OpenDomainStrategy.AS_IS)
#print (OpenDomainStrategy("0").name)
print("AS_IS" in OpenDomainStrategy)