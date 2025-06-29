import requests
from bs4 import BeautifulSoup
import json

def getHTML(url):
    result = requests.get(url)
    return result.text

def main():
    html = getHTML('https://it-events.com/')
    soup_html=BeautifulSoup(html, 'lxml')
    events_list=soup_html.findAll('a', {'class' : 'event-list-item__title'})
    #events_dates = soup_html.findAll('div', {'class': 'event-list-item__info'})

    for event in events_list:
        print(event.previousSibling.previousSibling)
    #print(events_list)
    print(requests.get('https://it-events.com/'))
    for event in events_list:
        print(event.text.strip(), '(', event.nextSibling.nextSibling.text.strip(), ')')

   # with open("events_list.json", "w") as file:
    #    json.dump(events_list, file, indent=4, ensure_ascii=False)
   # with open("events_list.json") as file:
    #    events_list=json.load(file)
    #print(events_list)

main()
#if __name__ =='__main__':       #если запускается как скрипт, то работает, если импортируется как модуль - ничего не делает
 #   main()