#define TIPO long int
#include <iostream>
#include <math.h>
#include <time.h>

struct node
{
  TIPO data;
  node *next;
};


class linked_list
{
  public:
    node *first, *last;
      int num;

  public:
    linked_list()
      { first = NULL;
        last = NULL;
         num = 0;
    }
    void insert(TIPO x)
      {
        node *p = new node;
      p->data = x;
        p->next = NULL;
        if(first == NULL)
        { first = p;
            last = p;
            num = 1;
         }
        else
        { last->next = p;
          last = p;
            num++;
          }
      }
      void free_linked_list()
      {  if(num > 0)
        { node *p1, *p2;
          p1 = first;
          while (p1 != NULL )
          { p2 = p1;
              p1 = p1->next;
              delete(p2);
          }
            first = NULL;
          last = NULL;
            num = 0;
         }
      }
      void get_linked_list(linked_list L)
      {  if(num > 0)
          free_linked_list();
        first = L.first;
        last = L.last;
         num = L.num;
      }
      void clone_linked_list(linked_list lis)
      {
        if(num > 0)
          free_linked_list();
      if(lis.first != NULL)
      {  node *p = lis.first;
          while(p != NULL)
          { insert(p->data);
              p = p->next;
            }
      }
      }
};
