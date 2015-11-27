#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char **argv)
{

	FILE *fi, *fo, *cmd;
	ssize_t read;
	char buff[800];
	char *line, *token, *next_token;
	int total_lines;
	int i = 0;
	size_t len = 0;

	fi = fopen("9x1.adj","r");
	fo = fopen("9x1vec.bin","wb");

	cmd = popen("wc -l 9x1.adj", "r");
	fgets(buff, sizeof(buff), cmd);
	total_lines = atoi(buff);

	while(fgets(buff, sizeof buff, fi)!=NULL)
	{
		i = 0;
		int num;
		token = strtok(buff, " ");
		while(i < total_lines)
		{
			num = atoi(token);
			if(i == num)
			{
				num = 1;
				token = strtok(NULL, " ");
			}
			else
			{
				num = 0;
			}
			char snum[5];
			snprintf(snum,5,"%d",num);
			printf("Writing %d to file\n",num);
			fwrite(&num,sizeof num, 1, fo);
			i+=1;
		}
	}

	fclose(fo);
	fclose(fi);
}
