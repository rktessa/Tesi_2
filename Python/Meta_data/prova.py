import requests
import json
import pandas

# request data from api
response = requests.get("https://covidmap.umd.edu/api/resources?indicator=covid&type=smoothed&country=France&daterange=20201115-20201130").text

#convert json data to dic data for use!
jsonData = json.loads(response)

# convert to pandas dataframe
df = pandas.DataFrame.from_dict(jsonData['data'])


# Get the list of all column names from headers

column_names = list(df.columns.values)

print(column_names)


# Other request tentative


response2 = requests.get("https://covidmap.umd.edu/apiv2/resources?indicator=covid&type=smoothed&country=Italy")
#convert json data to dic data for use!
jsonData2 = json.loads(response)

# convert to pandas dataframe
df2 = pandas.DataFrame.from_dict(jsonData['data'])


# Get the list of all column names from headers

column_names2 = list(df.columns.values)

print(column_names2)