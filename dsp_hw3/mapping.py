import sys

def read_file(file_name):
	mapping = dict()
	with open(file_name, 'r', encoding='big5-hkscs') as file:
		content = file.readlines()
		for line in content:
			line = line.split(' ')
			mapping[line[0]] = line[0]
			phone = line[1].split('/')
			for s in phone:
				if mapping.__contains__(s[0]):
					mapping[s[0]] += " "+line[0]
				else:
					mapping[s[0]] = line[0]

	return mapping
		
def write_file(file_name, dictionary):
	with open(file_name, 'w', encoding='big5-hkscs') as file:
		for key,value in dictionary.items():
			file.write('%s %s\n' % (key, value))

input_file = sys.argv[1]
target_file = sys.argv[2]
mapping = read_file(input_file)
write_file(target_file, mapping)
